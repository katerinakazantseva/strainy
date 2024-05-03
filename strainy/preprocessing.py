import subprocess
import os
import logging
import pysam
import gfapy

from strainy.params import StRainyArgs

logger = logging.getLogger()

def create_bam_file(fasta_file, fastq_file, output_file, num_threads, index=True):
    """
    Create a .bam file, requires user provided --fastq argument containing reads.
    """
    if not os.path.isfile(fastq_file):
        raise Exception("Reads file not found")

    logger.info(f"Creating bam file from {fasta_file} and {fastq_file}")
    minimap_mode = "map-ont" if StRainyArgs().mode == "nano" else "map-hifi"
    subprocess.check_output(f"minimap2 -ax {minimap_mode} {fasta_file} {fastq_file} -t {num_threads} | " \
                            f"samtools sort -@4 -t {num_threads} > {output_file}",
                            shell=True)
    if index:
        pysam.samtools.index(f"{output_file}", f"{output_file}.bai")
    logger.info(".bam file created!")


def gfa_to_fasta(gfa_file, output_file):
    """
    Creates a fasta file from the input gfa file. This is needed if the user
    omitted the optional -f argument or if the input graph is modified as a
    result of --split-long-unitigs argument.
    """
    fasta_cmd = f"""awk '/^S/{{print ">"$2"\\n"$3}}' {gfa_file} > {output_file}"""
    try:
        logger.debug(f"Creating fasta file from the gfa file {gfa_file}")
        subprocess.check_output(fasta_cmd, shell=True, capture_output=False, stderr=open(os.devnull, "w"))
        #logger.info("Done!")

    except subprocess.CalledProcessError as e:
        logger.error(e)
        raise Exception(f'Error while creating fasta file from the gfa file: {gfa_file} ' \
                        f'Optionally, you can create a fasta file yourself and provide it with "-f file.fasta"')


def add_gfa_line(input_graph, *args):
    """
    Add a gfa line to the input graph. Works for S and L lines.
    """
    args = [a for a in args if a is not None]
    line = "\t".join(args)
    try:
        input_graph.add_line(line)
    except gfapy.NotUniqueError:
        logger.warning(f"Tried insterting duplicate line, ignoring:\n{line}")
        pass

def split_long_unitigs(input_graph, output_file):
    """
    Replaces unitigs with lengths greater than StRainyArgs().splen with multiple
    shorter unitigs for better load balancing accross threads. The number of new
    unitigs will be created is the ceiling of unitig.length / StRainyArgs().splen
    The leftmost newly created unitig inherits the incoming (dovetails_L) edges
    of the original unitig and the rightmost new unitig inherits the outgoing
    (dovetails_R) edges. Other unitigs in between form a chain from one end to
    the other.
    The modified graph is saved to a new file and used for the rest of the
    stRainy pipeline inestead of the original graph. New bam and fasta files
    need to be generated.
    """
    # Convert user provided -splen from kb to b
    split_length = int(StRainyArgs().splen * 1000)
    for unitig in input_graph.segments:
        if unitig.length > split_length:
            # number of unitigs that will replace the original unitig
            n_new_unitigs = -(unitig.length // -split_length)
            # length of the each new unitig
            new_unitig_len = unitig.length // n_new_unitigs

            # edges of the original unitig that will be removed
            to_remove = []
            for i in range(n_new_unitigs):
                new_unitig_name = f"{unitig.name}_s{i+1}"
                if i == n_new_unitigs - 1:
                    # the last of the new unitigs gets the remaining bases too
                    new_unitig_seq = unitig.sequence[i*new_unitig_len:]
                else:
                    # all other unitigs get new_unitig_len number of bases
                    new_unitig_seq = unitig.sequence[i*new_unitig_len : (i+1) * new_unitig_len]

                # unitig string has the form ['S', name, sequence, ..fields..]
                optional_fields = str(unitig).split("\t")[3:]
                add_gfa_line(input_graph, "S", new_unitig_name, new_unitig_seq, *optional_fields)

                # leftmost new unitigs inherits the L edges
                if i == 0:
                    for edge in unitig.dovetails_L:
                        # new_edge_str[1]: from_name
                        # new_edge_str[3]: to_name
                        new_edge_str = str(edge).split("\t")
                        if edge.from_name == unitig.name:
                            new_edge_str[1] = new_unitig_name
                        else:
                            new_edge_str[3] = new_unitig_name
                        add_gfa_line(input_graph, *new_edge_str)
                        # store the original edge to be removed later on
                        to_remove.append(edge)
                else:
                    # connect the new unitigs to one another
                    add_gfa_line(input_graph,
                                 "L",
                                 prev_unitig, "+",
                                 new_unitig_name, "+",
                                 "0M")
                # The rightmost unitig inherits the R edges
                if i == n_new_unitigs - 1:
                    for edge in unitig.dovetails_R:
                        new_edge_str = str(edge).split("\t")
                        if edge.from_name == unitig.name:
                            new_edge_str[1] = new_unitig_name
                        else:
                            new_edge_str[3] = new_unitig_name
                        add_gfa_line(input_graph, *new_edge_str)
                        to_remove.append(edge)

                prev_unitig = new_unitig_name
            # remove the original unitig and its connections from the graph
            for e in to_remove:
                try:
                    input_graph.rm(e)
                except gfapy.error.RuntimeError:   #in case of self-loops
                    pass
            unitig.disconnect()

    for path in input_graph.paths:
        input_graph.rm(path)

    input_graph.to_file(output_file)


def get_unitigs_to_phase(input_graph, bam_file):
    """
    Returns a list of unitig names that are fit to phase based on the user defined
    min_unitig_coverage, max_unitig_coverage, min_unitig length parameters.
    """
    edges_to_phase = []
    min_unitig_length = 1000 * StRainyArgs().min_unitig_length # convert kb to b
    for unitig in input_graph.segments:
        alignment_coverage = round(float(pysam.samtools.coverage("-r",
                                                                 unitig.name,
                                                                 bam_file,
                                                                 "--no-header").
                                                                 split()[6]))
        if (StRainyArgs().min_unitig_coverage <= alignment_coverage <= StRainyArgs().max_unitig_coverage
                and unitig.length > min_unitig_length):
            edges_to_phase.append(unitig.name)

    return edges_to_phase


def preprocess_cmd_args(args):
    """
    Do preprocessing based on the input cmd arguments before starting phasing
    or transforming. Accessing arguments via args.XX instead of stRainyArguments.XX
    as some arguments may not be initialized yet.
    """

    if (args.bam or args.snp) and args.unitig_split_length != 0:
        logger.error("--bam and --snp arguments are incompatible with long edge splitting "
                     "(enabled by default). Add \"--unitig-split-length 0\" to disable..")
        raise Exception("Arguments exception")

    if args.snp and not args.bam:
        logger.error("--snp requires --bam to be set up")
        raise Exception("Arguments exception")

    if args.bam:
        if not os.path.isfile(args.bam):
            logger.error("bam file not found")
            raise Exception("Arguments exception")
        elif not os.path.isfile(args.bam + ".bai"):
            logger.error("bam file must be indexed (with samtools)")

    if args.snp:
        if not os.path.isfile(args.snp):
            logger.error("SNP/VCF file not found")
            raise Exception("Arguments exception")
        if not os.path.isfile(args.snp + ".tbi"):
            logger.error("SNP/VCF file must be block-gzipped and indexed (with bgzip + tabix)")
            raise Exception("Arguments exception")

    preprocessing_dir = os.path.join(args.output, "preprocessing_data")
    os.makedirs(preprocessing_dir, exist_ok=True)
    input_graph = gfapy.Gfa.from_file(args.gfa)
    if args.unitig_split_length != 0:
        split_long_unitigs(input_graph,
                           os.path.join(preprocessing_dir, "long_unitigs_split.gfa"))
        args.gfa = os.path.join(preprocessing_dir, "long_unitigs_split.gfa")
        args.graph_edges = input_graph.segment_names

    if args.fasta is None or args.unitig_split_length != 0:
        gfa_to_fasta(args.gfa,
                     os.path.join(preprocessing_dir,"gfa_converted.fasta"))
        args.fasta = os.path.join(preprocessing_dir,"gfa_converted.fasta")

    if args.bam is None or args.unitig_split_length != 0:
        create_bam_file(args.fasta,
                        args.fastq,
                        os.path.join(preprocessing_dir, "long_unitigs_split.bam"),
                        args.threads)
        args.bam = os.path.join(preprocessing_dir, "long_unitigs_split.bam")

    logger.info("Checking which sequences need to be phased")
    args.edges_to_phase = get_unitigs_to_phase(input_graph, args.bam)
    filtered_out = set(args.graph_edges) - set(args.edges_to_phase)
    logger.info(f"{len(filtered_out)}/{len(args.graph_edges)} unitigs will NOT be phased.")
    # args.graph_edges = args.edges_to_phase

    # log filtered out files
    with open(f'{StRainyArgs().output_intermediate}/filtered_out_unitigs.txt', 'w+') as f:
        for i in filtered_out:
            f.write(f'{i}\n')


#!/usr/bin/env python

import sys
import os
import subprocess
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import gfapy
import logging
import shutil

import pysam

from strainy.phase import phase_main
from strainy.transform import transform_main
from strainy.params import StRainyArgs, init_global_args_storage
from strainy.logging import set_thread_logging


logger = logging.getLogger()

def create_fasta_file(gfa_file):
    """
    Creates a fasta file from the input gfa file. This is needed if the user
    omitted the optional -f argument or if the input graph is modified as a
    result of --split-long-unitigs argument.
    """
    fasta_name = os.path.join(StRainyArgs().output, 'gfa_converted.fasta')
    fasta_cmd = f"""awk '/^S/{{print ">"$2"\\n"$3}}' {gfa_file} > {fasta_name}"""
    try:
        logger.info(f'Creating fasta file from the gfa file {gfa_file}')
        subprocess.check_output(fasta_cmd, shell=True, capture_output=False, stderr=open(os.devnull, "w"))
        StRainyArgs().fasta = fasta_name
        logger.info('Done!')

    except subprocess.CalledProcessError as e:
        logger.error(e)
        raise Exception(f'Error while creating fasta file from the gfa file: {gfa_file} ' \
                        f'Optionally, you can create a fasta file yourself and provide it with "-f file.fasta"')


def split_long_unitigs(input_graph):
    """
    Replaces unitigs with lengths greater than StRainyArgs().splen with multiple
    shorter unitigs for better load balancing accross threads. The number of new
    unitigs will be created is the ceiling of 
    unitig.length / StRainyArgs().splen
    The leftmost newly created unitig inherits the incoming (dovetails_L) edges
    of the original unitig and the rightmost new unitig inherits the outgoing
    (dovetails_R) edges. Other unitigs in between form a chain from one end to
    the other.
    The modified graph is saved to a new file and used for the rest of the
    stRainy pipeline inestead of the original graph. If there is at least one
    unitig in the original graph that is split, new bam and fasta files need to
    be generated. Fasta file is generated after exiting this function.
    """
    graph_modified = False
    split_length = int(StRainyArgs().splen * 1000) # convert kb to b
    for unitig in input_graph.segments:
        if unitig.length > split_length:
            graph_modified = True
            n_new_unitigs = -(unitig.length // -split_length) # ceiling
            new_unitig_len = unitig.length // n_new_unitigs
            for i in range(n_new_unitigs):
                new_unitig_name = f"{unitig.name}_s{i+1}"
                new_unitig_seq = unitig.sequence[i*new_unitig_len : (i+1) * new_unitig_len]
                input_graph.add_line(f"S\t{new_unitig_name}\t{new_unitig_seq}")
                
                # the left most unitig and the right most unitig
                # inherits the incoming edges of the original unitig
                if i == 0:
                    # the left most unitig inherits the incoming edges
                    for edge in unitig.dovetails_L:
                        if edge.from_name == unitig.name:
                            input_graph.add_line(
                                f"L"
                                f"\t{new_unitig_name}\t"
                                f"{edge.from_orient}\t"
                                f"{edge.to_name}\t"
                                f"{edge.to_orient}\t"
                                f"0M"
                                )
                        else:
                            input_graph.add_line(
                                f"L\t"
                                f"{edge.from_name}\t"
                                f"{edge.from_orient}\t"
                                f"{new_unitig_name}\t"
                                f"{edge.to_orient}\t"
                                f"0M"
                                )
                else:
                    # connect all the unitigs (except for the first one) to the previous unitig
                    input_graph.add_line(
                                f"L\t"
                                f"{prev_unitig}\t"
                                f"+\t"
                                f"{new_unitig_name}\t"
                                f"+\t"
                                f"0M"
                                )
                if i == n_new_unitigs - 1:
                    for edge in unitig.dovetails_R:
                        if edge.from_name == unitig.name:
                            input_graph.add_line(
                                f"L"
                                f"\t{new_unitig_name}\t"
                                f"{edge.from_orient}\t"
                                f"{edge.to_name}\t"
                                f"{edge.to_orient}\t"
                                f"0M"
                                )
                        else:
                            input_graph.add_line(
                                f"L\t"
                                f"{edge.from_name}\t"
                                f"{edge.from_orient}\t"
                                f"{new_unitig_name}\t"
                                f"{edge.to_orient}\t"
                                f"0M"
                                )
                prev_unitig = new_unitig_name

            unitig.disconnect()

    if graph_modified:
        # save the new graph
        modified_gfa_path = os.path.join(StRainyArgs().output,  
                                'long_unitigs_split.gfa')
        # new bam file needs to be created
        output_bam_path = os.path.join(StRainyArgs().output, 
                                'long_unitigs_split.bam')
        input_graph.to_file(modified_gfa_path)
        minimap_mode = "map-ont" if StRainyArgs().mode == "nano" else "map-hifi"
        # subprocess.check_output(f"minimap2 -ax {minimap_mode} {modified_gfa_path} {StRainyArgs().fastq} | " \
        #                         f"samtools sort -@4 -t {StRainyArgs.threads - 1} > {output_bam_path}",
        #                         shell=True)
        # pysam.samtools.index(f"{output_bam_path}", f"{output_bam_path}.bai")
        logger.info(f"--split-long-unitigs flag is set, modified gfa file ({modified_gfa_path}) " \
                    f"will be used for the computation")
        return True, modified_gfa_path, output_bam_path

    return False, '', ''


def main():
    #Setting executable paths
    strainy_root = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, strainy_root)

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", help="output dir",required=True)
    requiredNamed.add_argument("-b", "--bam", help="bam file",required=True)
    requiredNamed.add_argument("-g", "--gfa", help="gfa file",required=True)
    requiredNamed.add_argument("-m", "--mode", help="", choices=["hifi", "nano"], required=True)
    
    parser.add_argument("stage", help="stage to run: either phase or transform", choices=["phase", "transform"])
    parser.add_argument("-s", "--snp", help="vcf file", default=None)
    parser.add_argument("-t", "--threads", help="number of threads to use", type=int, default=4)
    parser.add_argument("-f", "--fasta", help="fasta file", required=False)
    # TODO: update the README.md file to include the new arguments
    parser.add_argument("-fq", "--fastq",
                        help="fastq file containing reads to perform alignment, only used if -splu is set",
                        required=False,
                        default=False)
    parser.add_argument("-splen", "--unitig-split-length",
                        help="The length (in kb) which the  unitigs that are longer will be split, only used if -splu is set",
                        required=False,
                        type=float,
                        default=50)
    parser.add_argument("-splu", "--split-long-unitigs",
                        help="split unitigs with long sequences into multiple smaller unitigs for faster processing, requires -fq to be provided",
                        action='store_true',
                        required=False)

    args = parser.parse_args()
    args.strainy_root = strainy_root

    bam_index = args.bam + ".bai"
    bam_index_exist = os.path.exists(bam_index)
    if bam_index_exist == False:
        raise Exception("No index file found (%s) Please create index using \"samtools index\"." % bam_index)

    #important so that global variables are inherited
    #multiprocessing.set_start_method("fork")

    #list of edges to process (can be manually modified for debugging)
    input_graph = gfapy.Gfa.from_file(args.gfa)
    args.graph_edges = input_graph.segment_names

    #setting up global arguments storage
    init_global_args_storage(args)

    if bool(StRainyArgs().fastq) != StRainyArgs().splu:
        parser.error("To split long unitigs, --fastq should be provided together" \
                        "with the --split-long-unitigs flag. Alternatively, you can omit both" \
                        "to run stRainy without without splitting the long unitigs.")

    BIN_TOOLS = ["samtools", "bcftools", StRainyArgs().flye]
    for tool in BIN_TOOLS:
        if not shutil.which(tool):
            print("{} not installed".format(tool), file=sys.stderr)
            return 1

    if not os.path.isdir(StRainyArgs().output):
        os.mkdir(StRainyArgs().output)

    set_thread_logging(StRainyArgs().output, "root", None)

    if StRainyArgs().splu:
        graph_modified, new_gfa_path, new_bam_path = split_long_unitigs(input_graph)
        if graph_modified:
            # If the graph is modified, need to update some stRainy arguments
            StRainyArgs().bam = new_bam_path
            StRainyArgs().gfa = new_gfa_path
            StRainyArgs().edges = input_graph.segment_names

    if graph_modified or args.fasta is None:
        create_fasta_file(StRainyArgs().gfa)

    if args.stage == "phase":
        sys.exit(phase_main(args))
    elif args.stage == "transform":
        sys.exit(transform_main(args))


if __name__ == "__main__":
    main()


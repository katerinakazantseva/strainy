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

def split_long_unitigs(input_graph, N=15000):
    # TODO: change N
    # splits unitigs S with a sequence longer than N into two new unitigs (S1, S2) connected to each other
    # incoming edges of S go to S1, and the outgoing edges goes from S2
    graph_modified = False
    for unitig in input_graph.segments:
        if unitig.length > N:
            graph_modified = True
            # prepare the data of the new unitigs
            # s1 takes the first half of the initial unitigs sequence, s2 takes the second half
            s1_sequence = unitig.sequence[: unitig.length // 2]
            s1_name = unitig.name + "_s1"
            s2_sequence = unitig.sequence[unitig.length // 2 :]
            s2_name = unitig.name + "_s2"

            # add the segments
            input_graph.add_line(f"S\t{s1_name}\t{s1_sequence}")
            input_graph.add_line(f"S\t{s2_name}\t{s2_sequence}")

            # link s1 to s2
            input_graph.add_line(f"L\t{s1_name}\t+\t{s2_name}\t+\t0M")

            # add incoming links
            for edge in unitig.dovetails_L:
                input_graph.add_line(f"L\t{edge.from_segment.name}\t+\t{s1_name}\t+\t0M")

            # add outgoing links
            for edge in unitig.dovetails_R:
                input_graph.add_line(f"L\t{s2_name}\t+\t{edge.to_segment.name}\t+\t0M")

            # finally, remove the original unitig that is replaced by s1 and s2
            unitig.disconnect()

    if graph_modified:
        # if the graph is modified, new bam file needs to be created
        modified_gfa_path = os.path.join(StRainyArgs().output,  
                                'long_unitigs_split.gfa')
        output_bam_path = os.path.join(StRainyArgs().output, 
                                'long_unitigs_split.bam')
        input_graph.to_file(modified_gfa_path)
        # TODO: make sure the minimap modes are correct and check the minimap command
        minimap_mode = "map-ont" if StRainyArgs().mode == "nano" else "map-hifi"
        # subprocess.check_output(f"minimap2 -ax {minimap_mode} {modified_gfa_path} {StRainyArgs().fastq} | " \
        #                         f"samtools sort -@4 -t {StRainyArgs.threads - 1} > {output_bam_path}",
        #                         shell=True)
        pysam.samtools.index(f"{output_bam_path}", f"{output_bam_path}.bai")
        logger.info(f"--split-long-unitigs flag is set, a new gfa file ({modified_gfa_path}) " \
                        f"which is created based on the input gfa file ({StRainyArgs().gfa}) will be used for the computation")
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
    parser.add_argument("-fq", "--fastq", help="fastq file containing reads to perform alignment, only used if -slu is set", required=False, default=False)
    parser.add_argument("-slu", "--split-long-unitigs",
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

    if bool(StRainyArgs().fastq) != StRainyArgs().slu:
        parser.error("To split long unitigs, --fastq file should be provided together" \
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

    if StRainyArgs().slu:
        modified, new_gfa_path, new_bam_path = split_long_unitigs(input_graph)
        if modified:
            # If the graph is modified, need to update some stRainy arguments
            StRainyArgs().bam = new_bam_path
            StRainyArgs().gfa = new_gfa_path
            StRainyArgs().edges = input_graph.segment_names
            # TODO: do I need to recalculate the fasta file?

    if args.fasta is None:
        fasta_name = os.path.join(StRainyArgs().output, 'gfa_converted.fasta')
        fasta_cmd = f"""awk '/^S/{{print ">"$2"\\n"$3}}' {StRainyArgs().gfa} > {fasta_name}"""
        try:
            logger.info(f'Creating fasta file from the gfa file {StRainyArgs().gfa}')
            subprocess.check_output(fasta_cmd, shell=True, capture_output=False, stderr=open(os.devnull, "w"))
            StRainyArgs().fasta = fasta_name
            logger.info('Done!')

        except subprocess.CalledProcessError as e:
            logger.error(e)
            logger.error(f'Error while creating fasta file from the gfa file: {StRainyArgs().gfa} '
                         f'Optionally, you can create a fasta file yourself and provide it with "-f file.fasta"')
            return 1

    if args.stage == "phase":
        sys.exit(phase_main(args))
    elif args.stage == "transform":
        sys.exit(transform_main(args))


if __name__ == "__main__":
    main()


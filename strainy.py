#!/usr/bin/env python

import sys
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
import shutil

import gfapy

from strainy.phase import phase_main
from strainy.transform import transform_main
from strainy.params import StRainyArgs, init_global_args_storage
from strainy.logging import set_thread_logging
from strainy.preprocessing import preprocess_cmd_args


logger = logging.getLogger()


def main():
    #Setting executable paths
    strainy_root = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, strainy_root)

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", help="output dir",required=True)
    requiredNamed.add_argument("-b", "--bam", help="bam file",required=True)
    requiredNamed.add_argument("-g", "--gfa", help="gfa file",required=True)
    requiredNamed.add_argument("-m", "--mode", help="", choices=["hifi", "nano"], required=True)
    
    parser.add_argument("stage", help="stage to run: either phase or transform", choices=["phase", "transform"])
    parser.add_argument("-s", "--snp", help="vcf file", default=None)
    parser.add_argument("-t", "--threads", help="number of threads to use", type=int, default=4)
    parser.add_argument("-f", "--fasta", help="fasta file", required=False)
    parser.add_argument("-splu", "--split-long-unitigs",
                        help="Split unitigs with long sequences into multiple smaller unitigs for faster processing, requires -fq to be provided",
                        action='store_true',
                        required=False)
    parser.add_argument("-fq", "--fastq",
                        help="fastq file containing reads to perform alignment, only used if -splu is set",
                        required=False,
                        default=False)
    parser.add_argument("-splen", "--unitig-split-length",
                        help="The length (in kb) which the unitigs that are longer will be split, only used if -splu is set",
                        required=False,
                        type=float,
                        default=50)


    args = parser.parse_args()
    args.strainy_root = strainy_root

    bam_index = args.bam + ".bai"
    bam_index_exist = os.path.exists(bam_index)
    if bam_index_exist == False:
        raise Exception("No index file found (%s) Please create index using \"samtools index\"." % bam_index)

    #important so that global variables are inherited
    #multiprocessing.set_start_method("fork")

    #setting up global arguments storage
    input_graph = gfapy.Gfa.from_file(args.gfa)
    args.graph_edges = input_graph.segment_names
    init_global_args_storage(args)

    BIN_TOOLS = ["samtools", "bcftools", StRainyArgs().flye]
    for tool in BIN_TOOLS:
        if not shutil.which(tool):
            print("{} not installed".format(tool), file=sys.stderr)
            return 1
        
    if not os.path.isdir(StRainyArgs().output):
        os.mkdir(StRainyArgs().output)
    set_thread_logging(StRainyArgs().output, "root", None)

    preprocess_cmd_args(args, parser)

    # set one more time for the modified args
    init_global_args_storage(args)
    
    if args.stage == "phase":
        sys.exit(phase_main(args))
    elif args.stage == "transform":
        sys.exit(transform_main(args))


if __name__ == "__main__":
    main()


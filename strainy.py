#!/usr/bin/env python

import sys
import os
import re
import subprocess
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import gfapy
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

    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument("-o", "--output", help="directory that will contain the output files",required=True)
    requiredNamed.add_argument("-g", "--gfa", help="gfa file",required=True)
    requiredNamed.add_argument("-m", "--mode", help="type of reads", choices=["hifi", "nano"], required=True)
    requiredNamed.add_argument("-q", "--fastq",
                    help="fastq file containing reads to perform alignment, used to create a .bam file",
                    required=True)
    
    parser.add_argument("-s", "--stage", help="stage to run: either phase, transform or e2e (phase + transform)", choices=["phase", "transform", "e2e"], default="e2e")
    parser.add_argument("--snp", help="vcf file", default=None)
    parser.add_argument("-t", "--threads", help="number of threads to use", type=int, default=4)
    parser.add_argument("-f", "--fasta", help="fasta file", required=False)
    parser.add_argument("-b", "--bam", help="bam file",required=False)
    parser.add_argument("--link-simplify", required=False, action="store_true", default=False, dest="link_simplify",
                        help="Enable agressive graph simplification")
    parser.add_argument("--unitig-split-length",
                        help="The length (in kb) which the unitigs that are longer will be split, set 0 to disable",
                        required=False,
                        type=int,
                        default=50)
    parser.add_argument("--only_split",help="Do not run stRainy, only split long gfa unitigs", default='False', required=False)

    args = parser.parse_args()
    args.strainy_root = strainy_root

    #setting up global arguments storage
    input_graph = gfapy.Gfa.from_file(args.gfa)
    args.graph_edges = input_graph.segment_names
    init_global_args_storage(args)

    BIN_TOOLS = ["samtools", "bcftools", "minimap2", StRainyArgs().flye]
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
    if args.only_split=='True':
        sys.exit()
    elif args.stage == "phase":
        sys.exit(phase_main(args))
    elif args.stage == "transform":
        sys.exit(transform_main(args))
    elif args.stage == "e2e":
        phase_main(args)
        logger.info("Phase stage completed, starting transform now...")
        transform_main(args)
        logger.info("Transform stage completed, exiting...")


if __name__ == "__main__":
    main()


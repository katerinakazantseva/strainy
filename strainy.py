#!/usr/bin/env python

import sys
import os
import re
import subprocess
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter
import gfapy
import multiprocessing
import logging
import shutil

from strainy.phase import phase_main
from strainy.transform import transform_main
from strainy.params import StRainyArgs, init_global_args_storage
from strainy.logging import set_thread_logging
import strainy.params as params


logger = logging.getLogger()


def main():
    #Setting executable paths
    strainy_root = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, strainy_root)

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("stage", help="stage to run: either phase or transform")
    parser.add_argument("-s", "--snp", help="vcf file", default=None)
    parser.add_argument("-t", "--threads", help="number of threads", type=int, default=4)
    parser.add_argument("-f", "--fasta", help="fasta file", required=False)

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", help="output dir",required=True)
    requiredNamed.add_argument("-b", "--bam", help="bam file",required=True)
    requiredNamed.add_argument("-g", "--gfa", help="gfa file",required=True)
    requiredNamed.add_argument("-m", "--mode", help="", choices=["hifi", "nano"], required=True)

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

    BIN_TOOLS = ["samtools", "bcftools", StRainyArgs().flye]
    for tool in BIN_TOOLS:
        if not shutil.which(tool):
            print("{} not installed".format(tool), file=sys.stderr)
            return 1

    if not os.path.isdir(StRainyArgs().output):
        os.mkdir(StRainyArgs().output)

    set_thread_logging(StRainyArgs().output, "root", None)

    if args.fasta is None:
        fasta_name = os.path.join(StRainyArgs().output, 'gfa_converted.fasta')
        fasta_cmd = f"""awk '/^S/{{print ">"$2"\\n"$3}}' {StRainyArgs().gfa} > {fasta_name}"""
        try:
            logger.info(f'Creating fasta file from the provided gfa file {StRainyArgs().gfa}')
            subprocess.check_output(fasta_cmd, shell=True, capture_output=False, stderr=open(os.devnull, "w"))
            args.fasta = fasta_name
            logger.info('Done!')

        except subprocess.CalledProcessError as e:
            logger.error(e)
            logger.error('You can create a fasta file yourself and provide it with "-f file.fasta"')
            logger.error(f'Error creating fasta file from the provided gfa file: {args.gfa}'
                         f'Optionally, you can create a fasta file yourself and provide it with "-f file.fasta"')
            return 1

    #setting up again (to update fasta parameter)
    init_global_args_storage(args)

    if args.stage == "phase":
        sys.exit(phase_main(args))
    elif args.stage == "transform":
        sys.exit(transform_main(args))
    else:
        raise Exception("Stage should be either phase or transform!")


if __name__ == "__main__":
    main()


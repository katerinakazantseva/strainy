#!/usr/bin/env python3

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
from strainy.params import StRainyArgs
from strainy.logging import set_thread_logging
import strainy.params as params


logger = logging.getLogger()


def main():
    #Setting executable paths
    strainy_root = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, strainy_root)

    BIN_TOOLS = ["samtools", "bcftools", params.flye]
    for tool in BIN_TOOLS:
        if not shutil.which(tool):
            print("{} not installed".format(tool), file=sys.stderr)
            return 1

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

    bam_index = re.sub(".bam",".bam.bai", args.bam)
    bam_index_exist = os.path.exists(bam_index)
    if bam_index_exist == False:
        raise Exception("No index file found (%s) Please create index using \"samtools index\"." % bam_index)

    #important so that global variables are inherited
    multiprocessing.set_start_method("fork")

    #global arguments storage

    StRainyArgs.output = args.output
    StRainyArgs.bam = args.bam
    StRainyArgs.gfa = args.gfa
    StRainyArgs.mode = args.mode
    StRainyArgs.snp = args.snp
    StRainyArgs.threads = args.threads
    StRainyArgs.gfa_transformed = "%s/transformed_before_simplification.gfa" % args.output
    StRainyArgs.gfa_transformed1 =  "%s/transformed_after_simplification.gfa" % args.output
    StRainyArgs.gfa_transformed2 = "%s/transformed_after_simplification_merged.gfa" % args.output
    StRainyArgs.log_phase = os.path.join(args.output, "log_phase")
    StRainyArgs.log_transform = os.path.join(args.output, "log_transform")

    if not os.path.isdir(StRainyArgs.output):
        os.mkdir(StRainyArgs.output)

    fasta_name = os.path.join(StRainyArgs.output, 'gfa_converted.fasta')
    fasta_cmd = f"""awk '/^S/{{print ">"$2"\\n"$3}}' {StRainyArgs.gfa} | fold > {fasta_name}"""
    try:
        subprocess.check_output(fasta_cmd, shell=True, capture_output=False, stderr=open(os.devnull, "w"))
        StRainyArgs.fa = fasta_name
    except subprocess.CalledProcessError as e:
        print(e)
        logger.error(f'Error creating fasta file from the provided gfa file: {args.gfa}'
                     f'Optionally, you can create a fasta file yourself and provide it with "-f file.fasta"')
        return 1

    if not os.path.isdir(StRainyArgs.output):
        os.mkdir(StRainyArgs.output)

    input_graph = gfapy.Gfa.from_file(args.gfa)
    StRainyArgs.edges = input_graph.segment_names
    ###

    set_thread_logging(StRainyArgs.output, "root", None)

    if args.stage == "phase":
        sys.exit(phase_main())
    elif args.stage == "transform":
        sys.exit(transform_main())
    else:
        raise Exception("Stage should be aither phase or transform!")


if __name__ == "__main__":
    main()


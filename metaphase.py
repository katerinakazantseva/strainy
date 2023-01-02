#!/usr/bin/env python3

import sys
import os
import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter
import gfapy
import multiprocessing
import logging
import shutil

from metaphase.phase import phase_main
from metaphase.transform import transform_main
from metaphase.params import MetaPhaseArgs


logger = logging.getLogger()


def main():
    #Setting executable paths
    metaphase_root = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, metaphase_root)

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("stage", help="stage to run: either phase or transform")
    parser.add_argument("-s", "--snp", help="vcf file", default=None)
    parser.add_argument("-t", "--threads", help="number of threads", type=int, default=4)

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", help="output dir",required=True)
    requiredNamed.add_argument("-b", "--bam", help="bam file",required=True)
    requiredNamed.add_argument("-g", "--gfa", help="gfa file",required=True)
    requiredNamed.add_argument("-f", "--fa", help="fa file",required=True)

    args = parser.parse_args()

    bam_index = re.sub(".bam",".bam.bai", args.bam)
    bam_index_exist = os.path.exists(bam_index)
    if bam_index_exist == False:
        raise Exception("No index file found (%s) Please create index using \"samtools index\"." % bam_index)

    #important so that global variables are inherited
    multiprocessing.set_start_method("fork")

    #global arguments storage
    MetaPhaseArgs.output = args.output
    MetaPhaseArgs.bam = args.bam
    MetaPhaseArgs.gfa = args.gfa
    MetaPhaseArgs.fa = args.fa
    MetaPhaseArgs.snp = args.snp
    MetaPhaseArgs.threads = args.threads
    MetaPhaseArgs.gfa_transformed = "%s/transformed_before_simplification.gfa" % args.output
    MetaPhaseArgs.gfa_transformed1 =  "%s/transformed_after_simplification.gfa" % args.output
    MetaPhaseArgs.gfa_transformed2 = "%s/transformed_after_simplification_merged.gfa" % args.output
    MetaPhaseArgs.log_phase = os.path.join(args.output, "log_phase")
    MetaPhaseArgs.log_transform = os.path.join(args.output, "log_transform")

    input_graph = gfapy.Gfa.from_file(args.gfa)
    MetaPhaseArgs.edges = input_graph.segment_names
    ###

    #main_log = os.path.join(MetaPhaseArgs.log_dir, "root.log")
    #_enable_logging(main_log, debug=True)

    if args.stage == "phase":
        sys.exit(phase_main())
    elif args.stage == "transform":
        sys.exit(transform_main())
    else:
        raise Exception("Stage should be aither phase or transform!")


if __name__ == "__main__":
    main()


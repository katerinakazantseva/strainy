import sys
import platform
import re
import argparse
import gfapy
import logging
import shutil
import subprocess
import os
from strainy.__version__ import __version__
from strainy.params import StRainyArgs, init_global_args_storage
from strainy.params import *
import strainy.graph_operations.gfa_ops as gfa_ops
from strainy.all_steps import all_main


def _version():
    return __version__

logger = logging.getLogger()

def get_processor_name():
    if platform.system() == "Windows":
        return platform.processor()
    elif platform.system() == "Darwin":
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + '/usr/sbin'
        command ="sysctl -n machdep.cpu.brand_string"
        return subprocess.check_output(command).strip()
    elif platform.system() == "Linux":
        command = "cat /proc/cpuinfo"
        all_info = subprocess.check_output(command, shell=True).decode().strip()
        for line in all_info.split("\n"):
            if "model name" in line:
                return re.sub( ".*model name.*:", "", line,1)
    return ""



def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument("-o", "--output", help="output directory", required=True)
    requiredNamed.add_argument("-g", "--gfa_ref", help="input reference gfa to uncollapse", required=False)
    requiredNamed.add_argument("-m", "--mode", help="type of reads", choices=["hifi", "nano"], required=True)
    requiredNamed.add_argument("-q", "--fastq",
                               help="fastq file with reads to phase / assemble",
                               required=True)

    parser.add_argument("-s", "--stage", help="stage to run: either phase, transform or e2e (phase + transform)",
                        choices=["phase", "transform", "e2e"], default="e2e")
    parser.add_argument("--snp", help="path to vcf file with SNP calls to use", default=None)
    parser.add_argument("-t", "--threads", help="number of threads to use", type=int, default=4)
    parser.add_argument("-f", "--fasta_ref", required=False, help="input reference fasta to phase")
    parser.add_argument("-b", "--bam", help="path to indexed alignment in bam format", required=False)
    parser.add_argument("--link-simplify", required=False, action="store_true", default=False, dest="link_simplify",
                        help="Enable agressive graph simplification")
    parser.add_argument("--debug", required=False, action="store_true", default=False,
                        help="Generate extra output for debugging")
    parser.add_argument("--unitig-split-length",
                        help="The length (in kb) which the unitigs that are longer will be split, set 0 to disable",
                        required=False,
                        type=int,
                        default=50)
    parser.add_argument("--only-split", help="Do not run stRainy, only split long gfa unitigs", default='False',
                        required=False)
    parser.add_argument("-d", "--cluster-divergence", help="cluster divergence", type=float, default=0, required=False)
    parser.add_argument("-a", "--allele-frequency", help="Set allele frequency for internal caller only (pileup)",
                        type=float, default=0.2, required=False)
    parser.add_argument("--min-unitig-length",
                        help="The length (in kb) which the unitigs that are shorter will not be phased",
                        required=False,
                        type=float,
                        default=1)
    parser.add_argument("--min-unitig-coverage",
                        help="The minimum coverage threshold for phasing unitigs, unitigs with lower coverage will not be phased",
                        required=False,
                        type=int,
                        default=20)
    parser.add_argument("--max-unitig-coverage",
                        help="The maximum coverage threshold for phasing unitigs, unitigs with higher coverage will not be phased",
                        required=False,
                        type=int,
                        default=500)
    parser.add_argument("-v", "--version", action="version", version=_version())

    args = parser.parse_args()

    if args.gfa_ref is None and args.fasta_ref is None:
        logger.error("Reference GFA or FASTA must be provided")
        raise Exception("Arguments exception")
    elif args.gfa_ref is None:
        input_graph = gfa_ops.fa_to_gfa(args.fasta_ref)
    else:
        input_graph = gfapy.Gfa.from_file(args.gfa_ref)

    args.graph_edges = input_graph.segment_names
    args.edges_to_phase = []
    init_global_args_storage(args)
    BIN_TOOLS = ["samtools", "bcftools", "minimap2"]
    for tool in BIN_TOOLS:
        if not shutil.which(tool):
            print("{} not installed".format(tool), file=sys.stderr)
            return 1

    os.makedirs(StRainyArgs().output, exist_ok=True)
    os.makedirs(StRainyArgs().output_intermediate, exist_ok=True)
    if os.path.isdir(StRainyArgs().log_phase):
        shutil.rmtree(StRainyArgs().log_phase)
    os.makedirs(StRainyArgs().log_phase, exist_ok=True)

    preprocessing_dir = os.path.join(args.output, "preprocessing_data")
    os.makedirs(preprocessing_dir, exist_ok=True)
    if args.fasta_ref is None or args.unitig_split_length != 0:
        gfa_ops.gfa_to_fasta(args.gfa_ref,
                     os.path.join(preprocessing_dir,"gfa_converted.fasta"))
        args.fasta_ref = os.path.join(preprocessing_dir,"gfa_converted.fasta")


    init_global_args_storage(args)

    sys.exit(all_main(args))
    #if args.only_split=='True':
        #sys.exit()
    #elif args.stage == "phase":
    #    sys.exit(phase_main(args))
    #elif args.stage == "transform":
    #    sys.exit(transform_main(args))
    #elif args.stage == "e2e":







import pandas as pd
import gfapy

import pysam
import re
import os
import subprocess
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter

# Parse command line arguments
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-o", "--output", help="output dir",required = True)
requiredNamed.add_argument("-b", "--bam", help="bam file",required = True)
requiredNamed.add_argument("-g", "--gfa", help="gfa file",required = True)
requiredNamed.add_argument("-f", "--fa", help="fa file",required = True)
parser.add_argument("-s", "--snp", help="vcf file", default=None)
args = vars(parser.parse_args())

# Set up parameters

output = args["output"]
bam = args["bam"]
gfa = args["gfa"]
fa = args["fa"]
snp = args["snp"]


bam_index=re.sub(".bam",".bam.bai", bam)
bam_index_exist = os.path.exists(bam_index)
if bam_index_exist == False:
    raise Exception("No index file found (%s) Please create index using \"samtools index\"." % bam_index)

#transformed gfa file path
gfa_transformed = "%s/transformed_before_simplification.gfa" % output
gfa_transformed1 =  "%s/transformed_after_simplification.gfa" % output
gfa_transformed2 = "%s/transformed_after_simplification_merged.gfa" % output

minigraph=False

unseparated_cluster_min_reads=2

# Path to the installed Flye executable
# Should be Flye/bin/flye
flye = "../Flye/bin/flye"

# Path to store and read the consensus dictionary
# If one already exists and write_consensus_cache is true, it may be overwritten
consensus_cache_path = "consensus_dict_2.pkl"

# Whether to store the consensus cache or not
# This needs to be True to carry out it from phase part to transform part
write_consensus_cache = True

# Number processes, default (-1) uses all available
processes = -1

delete_files = False

# snp file path. If=None, metaPhase call snp using bcftools
snp = None

"""It is not recommended to change parameters below"""

# reads max mismatch count
R = 1


#cluster parameters
min_cluster_size=2
unclustered_group_N=1000000
unclustered_group_N2=3000000
split_id=10000

#transformation
parental_min_coverage=6
parental_min_len=0.7
start_end_gap=5
strong_cluster_min_reads=2
max_hops=10
min_reads_neighbour=0
min_reads_cluster=2

#simplification
cov_ratio=1.6



g = gfapy.Gfa.from_file(gfa)
edges=g.segment_names



# reads min intersection
I = 1000


# alignment filtering
clipp = 100
min_mapping_quality = 20
min_base_quality = 0
min_al_len = 1000
de_max = 0.05


# SNP allele frequency
AF = 0.1


# Please do not change
g = gfapy.Gfa.from_file(gfa)
edges = g.segment_names


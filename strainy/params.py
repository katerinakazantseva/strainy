import os


#TODO: this is a temporary way to share some global "constants"
#that are set at runtime once. Instead of adding one more constant,
#consider explicitly sharing it instead (or as a part of parsed command line arguments)
class GlobalArgumentStorage(object):
    pass

_glob_args = None

def init_global_args_storage(args):
    global _glob_args
    _glob_args = GlobalArgumentStorage()
    _glob_args.output = args.output
    _glob_args.bam = args.bam
    _glob_args.gfa = args.gfa
    _glob_args.mode = args.mode
    _glob_args.snp = args.snp
    _glob_args.threads = args.threads
    _glob_args.flye = os.path.join(args.strainy_root, "submodules", "Flye", "bin", "flye")
    #_glob_args.gfa_transformed = "%s/transformed_before_simplification.gfa" % args.output
    #_glob_args.gfa_transformed1 =  "%s/transformed_after_simplification.gfa" % args.output
    #_glob_args.gfa_transformed2 = "%s/transformed_after_simplification_merged.gfa" % args.output
    _glob_args.log_phase = os.path.join(args.output, "log_phase")
    _glob_args.log_transform = os.path.join(args.output, "log_transform")
    _glob_args.edges = args.graph_edges
    _glob_args.fa = args.fasta
    _glob_args.fq = args.fastq
    _glob_args.splen = args.unitig_split_length


def StRainyArgs():
    global _glob_args
    return _glob_args
##########


#TODO: constant storage

# Path to store and read the consensus dictionary
# If one already exists and write_consensus_cache is true, it may be overwritten
consensus_cache_path = "consensus_dict.pkl"

# Whether to store the consensus cache or not
# This needs to be True to carry out it from phase part to transform part
write_consensus_cache = True
delete_files = False

"""It is not recommended to change parameters below"""


#cluster parameters
R = 1 # reads max mismatch count
I = 1000 # reads min intersection
unseparated_cluster_min_reads = 3
min_cluster_size = 2
UNCLUSTERED_GROUP_N = 1000000
UNCLUSTERED_GROUP_N2 = 3000000
SPLIT_ID = 10000

#creating new unitigs
parental_min_coverage = 6
parental_min_len = 0.7
start_end_gap = 500
strong_cluster_min_reads = 2

#adding new links
max_hops = 10
min_reads_neighbour = 3
min_reads_cluster = 3

#simplification
cov_ratio = 1.6
minigraph = False

# alignment filtering
max_clipping = 100
min_mapping_quality = 20
min_base_quality = 0
min_al_len = 1000
# extended_aln_flank = 50
de_max = {"hifi": 0.05, "nano": 0.10}
min_consensus_cov = {"hifi": 3, "nano": 5}

# SNP allele frequency
AF = 0.1
split_allele_freq = 0.3

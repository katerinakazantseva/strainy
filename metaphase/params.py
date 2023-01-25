
class GlobalArgumentStorage(object):
    def __init__(self):
        self.output = None
        self.bam = None
        self.gfa = None
        self.fa = None
        self.snp = None
        self.edges = None

        self.gfa_transformed = None
        self.gfa_transformed1 =  None
        self.gfa_transformed2 = None

MetaPhaseArgs = GlobalArgumentStorage()

#TODO: constant storage

#TODO: link Flye as submodule

# Path to the installed Flye executable
# Should be Flye/bin/flye
flye = "/home/mkolmogo/projects/Flye/bin/flye"

minigraph=False
unseparated_cluster_min_reads=2

# Path to store and read the consensus dictionary
# If one already exists and write_consensus_cache is true, it may be overwritten
consensus_cache_path = "consensus_dict_2.pkl"

# Whether to store the consensus cache or not
# This needs to be True to carry out it from phase part to transform part
write_consensus_cache = True

# Number processes, default (-1) uses all available
#processes = 4

delete_files = False

# snp file path. If=None, metaPhase call snp using bcftools
# snp = None

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

#g = gfapy.Gfa.from_file(gfa)
#edges=g.segment_names

# reads min intersection
I = 1000

# alignment filtering
clipp = 100
min_mapping_quality = 20
min_base_quality = 0
min_al_len = 1000
extended_aln_flank = 5
de_max = {"hifi": 0.05, "nano": 0.10}

# SNP allele frequency
AF = 0.1

# Please do not change
#g = gfapy.Gfa.from_file(gfa)
#edges = g.segment_names

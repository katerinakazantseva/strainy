import pandas as pd
import gfapy

"""Please specify parameters below"""

# bam file path
bam = "/Users/ataberk/Documents/metagenomic-phasing/metaphase_data/results/sim3-3ecoli-alignment.bam"

# gfa file path
gfa = "/Users/ataberk/Documents/metagenomic-phasing/metaphase_data/flye_3ecoli_sim_noalt_haplo/assembly_graph.gfa"

# gaf file path
gaf_file = "/Users/ataberk/Documents/metagenomic-phasing/metaphase_data/results/flye_3ecoli_sim_noalt_haplo-hifi.sim.3ecoli.gaf"

# Where to save the transformed gfa file
gfa_transformed = "./3ecoli_modified_full.gfa"

# Path to the installed Flye executable
# Should be Flye/bin/flye
flye = "/Users/ataberk/Documents/metagenomic-phasing/software/Flye/bin/flye"

# Path to store and read the consensus dictionary
# If one already exists and write_consensus_cache is true, it may be overwritten
consensus_cache_path = "consensus_dict_2.pkl"

# Whether to store the consensus cache or not
# This needs to be True to carry out it from phase part to transform part
write_consensus_cache = True

# snp file path. If=None, metaPhase call snp using bcftools
snp = None

"""It is not recommended to change parameters below"""

# reads max mismatch count
R = 1

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
gf = pd.read_csv(gaf_file, sep="\t")
gf1 = pd.read_csv(gaf_file, sep=" ")
gaf = pd.concat((gf1[gf1.columns[0]], gf[gf.columns[5]],), axis=1, keys=['ReadName', 'Al'])
g = gfapy.Gfa.from_file(gfa)
edges = g.segment_names

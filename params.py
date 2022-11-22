import pandas as pd
import gfapy
import pysam


"""Please specify parameters below"""


#bam file path
bam=""

#gfa file path
gfa=""

#transformed gfa file path
gfa_transformed = ""


minigraph=False

#snp file path. If=None, metaPhase call snp using bcftools
snp=None


"""It is not recommended to change parameters below"""

#reads max mismatch count
R=1

#reads min intersection
I=1000

#alignment filtering
clipp=100
min_mapping_quality=20
min_base_quality=0
min_al_len=1000
de_max=0.05

#SNP allele frequency
AF=0.1

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


#Please do not change
gf=pd.read_csv(gaf_file, sep="\t")
gf1=pd.read_csv(gaf_file, sep=" ")
gaf= pd.concat((gf1[gf1.columns[0]],gf[gf.columns[5]],),axis=1, keys=['ReadName', 'Al'])
g = gfapy.Gfa.from_file(gfa)
edges=g.segment_names




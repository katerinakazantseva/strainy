import pandas as pd
import gfapy
import pysam


"""Please specify parameters below"""

#bam file path
bam="/Users/ekaterina.kazantseva/MT/5sim_staggered/flye_5ecoli_staggered_noalt_bubbles_005/5.bam"

#gfa file path
gfa = "/Users/ekaterina.kazantseva/MT/5sim_staggered/flye_5ecoli_staggered_noalt_bubbles_005/assembly_graph.gfa"

#gaf file path
gaf_file="/Users/ekaterina.kazantseva/Downloads/aln5.gaf"

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


#Please do not change
gf=pd.read_csv(gaf_file, sep="\t")
gf1=pd.read_csv(gaf_file, sep=" ")
gaf= pd.concat((gf1[gf1.columns[0]],gf[gf.columns[5]],),axis=1, keys=['ReadName', 'Al'])
g = gfapy.Gfa.from_file(gfa)
edges=g.segment_names
#edges=['edge_104', 'edge_310']
#edges=['edge_274', 'edge_387', 'edge_44','edge_293','edge_46']
#edges=['edge_387']



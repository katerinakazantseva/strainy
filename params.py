import pandas as pd
import gfapy
import pysam


"""Please specify parameters below"""

#bam file path
bam = "/Users/ataberk/Documents/metagenomic-phasing/metaphase_data/results/sim3-3ecoli-alignment.bam"

#gfa file path
gfa = "/Users/ataberk/Documents/metagenomic-phasing/metaphase_data/flye_3ecoli_sim_noalt_haplo/assembly_graph.gfa"

#gaf file path
gaf_file= "/Users/ataberk/Documents/metagenomic-phasing/metaphase_data/results/flye_3ecoli_sim_noalt_haplo-hifi.sim.3ecoli.gaf"

#transformed gfa file path
gfa_transformed = "./3ecoli_modified_3_edged_optimized.gfa"

#flye path
flye='/Users/ataberk/Documents/metagenomic-phasing/software/Flye/bin/flye'

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
edges=['edge_1', 'edge_2', 'edge_3']



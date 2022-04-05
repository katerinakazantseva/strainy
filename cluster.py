import csv
import pysam
from Bio import SeqIO
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mt
import sys,os,subprocess
from cluster_postprocess import postprocess
from build_adj_matrix import *
from community_detection import find_communities
#from color_bam import write_bam
#from color_bam import index

#from cdlib import algorithms
#bam="/Users/ekaterina.kazantseva/MT/test_data/SRR13128014_labeled.bam"
bam="/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/sim.3ecoliY.bam"
bam="/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/sim.3ecoliv2.bam"
#snp="/Users/ekaterina.kazantseva/MT/test.vcf"
snp=None

#edges=['edge_1', 'edge_10', 'edge_103', 'edge_105', 'edge_11', 'edge_111', 'edge_118', 'edge_124', 'edge_131', 'edge_143', 'edge_145', 'edge_146', 'edge_147', 'edge_149', 'edge_150', 'edge_151', 'edge_161', 'edge_162', 'edge_163', 'edge_174', 'edge_175', 'edge_176', 'edge_178', 'edge_181', 'edge_186', 'edge_187', 'edge_188', 'edge_19', 'edge_191', 'edge_192', 'edge_193', 'edge_194', 'edge_195', 'edge_196', 'edge_197', 'edge_198', 'edge_199', 'edge_2', 'edge_200', 'edge_201', 'edge_202', 'edge_203', 'edge_207', 'edge_208', 'edge_209', 'edge_214', 'edge_215', 'edge_220', 'edge_221', 'edge_222', 'edge_223', 'edge_230', 'edge_238', 'edge_24', 'edge_243', 'edge_244', 'edge_249', 'edge_25', 'edge_255', 'edge_258', 'edge_259', 'edge_26', 'edge_260', 'edge_261', 'edge_262', 'edge_263', 'edge_268', 'edge_269', 'edge_270', 'edge_271', 'edge_272', 'edge_277', 'edge_278', 'edge_279', 'edge_284', 'edge_285', 'edge_286', 'edge_287', 'edge_289', 'edge_292', 'edge_293', 'edge_294', 'edge_298', 'edge_3', 'edge_311', 'edge_312', 'edge_315', 'edge_316', 'edge_319', 'edge_320', 'edge_323', 'edge_325', 'edge_326', 'edge_327', 'edge_328', 'edge_329', 'edge_330', 'edge_331', 'edge_34', 'edge_343', 'edge_355', 'edge_36', 'edge_37', 'edge_372', 'edge_4', 'edge_40', 'edge_41', 'edge_42', 'edge_43', 'edge_44', 'edge_5', 'edge_543', 'edge_598', 'edge_599', 'edge_6', 'edge_600', 'edge_601', 'edge_615', 'edge_638', 'edge_656', 'edge_661', 'edge_662', 'edge_663', 'edge_667', 'edge_668', 'edge_669', 'edge_67', 'edge_670', 'edge_672', 'edge_673', 'edge_68', 'edge_680', 'edge_69', 'edge_694', 'edge_7', 'edge_70', 'edge_71', 'edge_72', 'edge_73', 'edge_74', 'edge_75', 'edge_765', 'edge_8', 'edge_85', 'edge_86', 'edge_89', 'edge_9', 'edge_90', 'edge_91']
edges=['edge_1', 'edge_10', 'edge_100', 'edge_101', 'edge_102', 'edge_103', 'edge_104', 'edge_105', 'edge_106', 'edge_107', 'edge_108', 'edge_109', 'edge_11', 'edge_110', 'edge_111', 'edge_112', 'edge_113', 'edge_114', 'edge_115', 'edge_116', 'edge_12', 'edge_13', 'edge_14', 'edge_15', 'edge_16', 'edge_17', 'edge_18', 'edge_19', 'edge_2', 'edge_20', 'edge_21', 'edge_22', 'edge_23', 'edge_24', 'edge_25', 'edge_26', 'edge_27', 'edge_28', 'edge_29', 'edge_3', 'edge_30', 'edge_31', 'edge_32', 'edge_33', 'edge_34', 'edge_35', 'edge_36', 'edge_37', 'edge_38', 'edge_39', 'edge_4', 'edge_40', 'edge_41', 'edge_42', 'edge_43', 'edge_44', 'edge_45', 'edge_46', 'edge_47', 'edge_48', 'edge_49', 'edge_5', 'edge_50', 'edge_51', 'edge_52', 'edge_53', 'edge_54', 'edge_55', 'edge_56', 'edge_57', 'edge_58', 'edge_59', 'edge_6', 'edge_60', 'edge_61', 'edge_62', 'edge_63', 'edge_64', 'edge_65', 'edge_66', 'edge_67', 'edge_68', 'edge_69', 'edge_7', 'edge_70', 'edge_71', 'edge_72', 'edge_73', 'edge_74', 'edge_75', 'edge_76', 'edge_77', 'edge_78', 'edge_79', 'edge_8', 'edge_80', 'edge_81', 'edge_82', 'edge_83', 'edge_84', 'edge_85', 'edge_86', 'edge_87', 'edge_88', 'edge_89', 'edge_9', 'edge_90', 'edge_91', 'edge_92', 'edge_93', 'edge_94', 'edge_95', 'edge_96', 'edge_97', 'edge_98', 'edge_99']

#edges=['edge_34']

R=1
I=1000
clipp=100
min_mapping_quality=20
min_base_quality=0
min_al_len=1000
de_max=0.05
AF=0.1



dirs = ("output/vcf/","output/adj_M/","output/clusters/","output/graphs/")
for dir in dirs:
    try:
        os.stat(dir)
    except:
         os.mkdir(dir)


def read_snp(snp):
    SNP_pos = []
    if snp==None:
        snpos = 'bcftools mpileup -r {} {} --no-reference -I --no-version --annotate FORMAT/AD | bcftools query -f  "%CHROM %POS [ %AD %DP]\n" >output/vcf/vcf_{}.txt'.format(edge,bam,edge)
        #snpos = 'bcftools mpileup -r {}  /Users/ekaterina.kazantseva/MT/test_data/SRR13128014_labeled.bam   --no-reference -I --no-version --annotate FORMAT/AD | bcftools query -f  "%CHROM %POS [ %AD %DP]\n" >output/vcf/vcf_{}.txt'.format(edge,edge)
        subprocess.check_output(snpos, shell=True, capture_output=False)
        with open("output/vcf/vcf_%s.txt" % edge) as f:
            lines = f.readlines()
            for line in lines:
                try:
                    AlFreq=int(str(line.split()[2]).split(',')[2])/int(line.split()[3])
                except(IndexError):  AlFreq=0
                if AlFreq>AF:
                    SNP_pos.append(line.split()[1])

    else:
        vcf = open(snp, "rt")
        for line in vcf:
            if line.split()[0] == edge:
                SNP_pos.append(line.split()[1])
    print(str(len(SNP_pos)) + " SNPs found")
    #print(SNP_pos)
    return(SNP_pos)


def read_bam(bam,SNP_pos, clipp=clipp):
    bamfile = pysam.AlignmentFile(bam, "rb")
    data = {}
    ln = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4]
    for pos in SNP_pos:
        for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools', min_base_quality=0,
                                           ignore_overlaps=False,min_mapping_quality=20,
                                           ignore_orphans=False, truncate=True):

            for pileupread in pileupcolumn.pileups:

                clipping=0
                start = pileupread.alignment.get_reference_positions()[0]
                stop = pileupread.alignment.get_reference_positions()[
                    len(pileupread.alignment.get_reference_positions()) - 1]

                de=float(str(pileupread).split('\t')[11].split('), (')[8].split(',')[1])
                for i in pileupread.alignment.cigartuples:
                    if i[0] == 4 or i[0] == 5:
                        if i[1]>clipp:
                            clipping=1


                if (clipping==1 or (stop-start)<min_al_len or de>de_max) and (int(start)!=0 and int(stop)!=int(ln)-1):
                    continue


                else:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        try:
                            data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[
                                pileupread.query_position]

                        except (KeyError):
                            data[pileupread.alignment.query_name] = {}
                            data[pileupread.alignment.query_name]["Start"] = pileupread.alignment.get_reference_positions()[
                            0]
                            data[pileupread.alignment.query_name]["Stop"] = pileupread.alignment.get_reference_positions()[
                                len(pileupread.alignment.get_reference_positions()) - 1]

                            data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[
                             pileupread.query_position]

    bamfile.close()
    return(data)


def clusters_vis_stats ( G,cl, clN,uncl, SNP_pos,bam, edge, I, AF,create_bam=False):
    cl.loc[cl['Cluster'] == 'NA', 'Cluster'] = 0
    cmap = plt.get_cmap('viridis')
    clusters = set(cl['Cluster'])
    cmap = cmap(np.linspace(0, 1, len(clusters)))
    colors = {}
    try:
        clusters.remove('0')
    except:
        KeyError
    colors[0] = "#505050"
    i = 0
    for cluster in clusters:
        colors[cluster] = mt.colors.to_hex(cmap[i])
        i = i + 1

    for index in cl.index:
        cl.loc[index, 'Color'] = colors[cl.loc[index, 'Cluster']]
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    nx.draw(G, nodelist=G.nodes(), with_labels=False, width=0.03, node_size=3, font_size=5,
            node_color=cl['Color'])
    # labels = cl['Cluster']

    ln = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4]
    cov = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[6]
    plt.suptitle(str(edge) + " coverage:" + str(cov) + " length:" + str(ln) + " clN:" + str(clN))
    plt.savefig("output/graphs/graph_%s_%s_%s.png" % (edge, I, AF), format="PNG", dpi=300)
    plt.close()
    # Calculate statistics
    print("Summary for: " + edge)
    print("Clusters found: " + str(clN))
    print("Reads unclassified: " + str(uncl))
    print("Number of reads in each cluster: ")
    print(cl['Cluster'].value_counts(dropna=False))

    stats = open('output/stats.txt', 'a')
    stats.write(edge + "\t" + str(ln) + "\t" + str(cov) + "\t" + str(len(cl['ReadName'])) + "\t" + str(
        len(SNP_pos)) + "\t" + str(clN) + "\t" + str(uncl) + "\n")
    stats.close()

    #if create_bam==True:
        #write_bam(bam, edge, I, AF)
        #index("output/bam/coloredBAM_%s_%s_%s.bam" % (edge, I, AF),edge, I, AF)


def main(edge):
    #READ SNPs
    print ("### Reading SNPs...")
    SNP_pos=read_snp(snp)
    if len(SNP_pos)==0:
        return

    #READ READS AND POSITIONS
    print ("### Reading Reads...")
    data=read_bam(bam,SNP_pos)
    cl=pd.DataFrame(data={'ReadName': data.keys()})
    cl['Cluster'] = 'NA'
    print (str(len(cl['ReadName']))+" reads found")
    if len(cl['ReadName'])==0:
        return

    #CALCULATE DISTANCE and ADJ MATRIX
    print ("### Calculatind distances/Building adj matrix...")
    #m=build_adj_matrix2(cl, data, SNP_pos, I, bam, edge, R)

    #m.to_csv("output/adj_M/adj_M_%s_%s_%s.csv" % (edge, I, AF))
    m=pd.read_csv("output/adj_M/adj_M_%s_%s_%s.csv" % (edge, I, AF),index_col='ReadName')
    print("### Removing overweighed egdes...")
    m=remove_edges (m, R)


    #BUILD graph and find clusters
    print("### Creating graph...")
    m1=m
    m1.columns = range(0,len(cl['ReadName']))
    m1.index=range(0,len(cl['ReadName']))
    G = nx.from_pandas_adjacency(change_w(m.transpose(), R))
    print("### Searching clusters...")
    cluster_membership = find_communities(G)
    clN = 0
    uncl = 0

    for value in set(cluster_membership.values()):
        group = [k for k, v in cluster_membership.items() if v == value]
        if len(group) > 3:
            clN = clN + 1
            cl['Cluster'][group] = value
        else:
            uncl = uncl + 1

    print(str(clN)+" clusters found")
    cl.to_csv("output/clusters/clusters_before_joining_%s_%s_%s.csv" % (edge, I, AF))

    print("### Cluster post-processing...")
    if clN!=0:
        cl=postprocess(bam, cl, SNP_pos, data, edge, R, I)
    clN=len(set(cl.loc[cl['Cluster']!='NA']['Cluster'].values))
    print(str(clN) + " clusters after post-processing")
    cl.to_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF))


    print("### Graph viz...")
    clusters_vis_stats (G,cl, clN,uncl,SNP_pos, bam, edge, I, AF,create_bam=False)



stats = open('output/stats.txt', 'a')
stats.write("Edge" + "\t" + "Len" + "\t" + "Coverage" + "\t" + "ReadsN" + "\t"+"SNPN"+"\t"+"ClustersN"+"\t"+"UnclusteredRN"+"\n")
stats.close()
for edge in edges:
    print("------------- "+str(edge)+" ------------- ")
    main(edge)



























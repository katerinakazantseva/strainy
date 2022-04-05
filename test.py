import csv
import pysam
from Bio import SeqIO
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from karateclub import LabelPropagation
from networkx.algorithms import community
from build_adj_matrix import build_adj_matrix2
from build_adj_matrix import *
from cluster_postprocess import cluster_consensuns
from cluster_postprocess import build_data
import sys,os,subprocess
import pygraphviz as gv
import re
#import pylab

edges=['edge_1', 'edge_10', 'edge_100', 'edge_101', 'edge_102', 'edge_103', 'edge_104', 'edge_105', 'edge_106', 'edge_107', 'edge_108', 'edge_109', 'edge_11', 'edge_110', 'edge_111', 'edge_112', 'edge_113', 'edge_114', 'edge_115', 'edge_116', 'edge_12', 'edge_13', 'edge_14', 'edge_15', 'edge_16', 'edge_17', 'edge_18', 'edge_19', 'edge_2', 'edge_20', 'edge_21', 'edge_22', 'edge_23', 'edge_24', 'edge_25', 'edge_26', 'edge_27', 'edge_28', 'edge_29', 'edge_3', 'edge_30', 'edge_31', 'edge_32', 'edge_33', 'edge_34', 'edge_35', 'edge_36', 'edge_37', 'edge_38', 'edge_39', 'edge_4', 'edge_40', 'edge_41', 'edge_42', 'edge_43', 'edge_44', 'edge_45', 'edge_46', 'edge_47', 'edge_48', 'edge_49', 'edge_5', 'edge_50', 'edge_51', 'edge_52', 'edge_53', 'edge_54', 'edge_55', 'edge_56', 'edge_57', 'edge_58', 'edge_59', 'edge_6', 'edge_60', 'edge_61', 'edge_62', 'edge_63', 'edge_64', 'edge_65', 'edge_66', 'edge_67', 'edge_68', 'edge_69', 'edge_7', 'edge_70', 'edge_71', 'edge_72', 'edge_73', 'edge_74', 'edge_75', 'edge_76', 'edge_77', 'edge_78', 'edge_79', 'edge_8', 'edge_80', 'edge_81', 'edge_82', 'edge_83', 'edge_84', 'edge_85', 'edge_86', 'edge_87', 'edge_88', 'edge_89', 'edge_9', 'edge_90', 'edge_91', 'edge_92', 'edge_93', 'edge_94', 'edge_95', 'edge_96', 'edge_97', 'edge_98', 'edge_99']
edges=[ 'edge_10', 'edge_101',  'edge_103', 'edge_104', 'edge_105', 'edge_108', 'edge_11','edge_110', 'edge_111',  'edge_114', 'edge_12', 'edge_15',  'edge_2', 'edge_21', 'edge_22', 'edge_23','edge_28', 'edge_29', 'edge_31',  'edge_35', 'edge_36','edge_4','edge_42', 'edge_47', 'edge_48',  'edge_5','edge_54','edge_60',  'edge_62',  'edge_65', 'edge_71', 'edge_72',  'edge_74', 'edge_75', 'edge_79', 'edge_84',  'edge_88', 'edge_89', 'edge_90', 'edge_91', 'edge_92', 'edge_94', 'edge_96', 'edge_97', 'edge_98', 'edge_99']

#edges=['edge_98','edge_101', 'edge_102' ]


#edges=['edge_98']
R=1
I=1000
clipp=100
min_mapping_quality=20
min_base_quality=0
min_al_len=1000
de_max=0.05
AF=0.1

bam="/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/sim.3ecoliv2.bam"

def read_snp(snp):
    SNP_pos = []
    if snp==None:
        snpos = 'bcftools mpileup -r {} {} --no-reference -I --no-version --annotate FORMAT/AD | bcftools query -f  "%CHROM %POS [ %AD %DP]\n" >output/vcf/vcf_{}.txt'.format(edge,bam,edge)
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
    #print(str(len(SNP_pos)) + " SNPs found")
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


def add_child_edge(edge,clN,g,cl,SNP_pos, data):
    f = g.segments
    #dat = {0: "A", 1: "A", 4: "A"}
    cons = build_data(cl, SNP_pos, data)
    cl_consensuns=cluster_consensuns(cl,cluster,SNP_pos, data, cons)
    for i in f:
        if i == edge:
            for key, val in cl_consensuns[clN].items():
                seq = i.sequence
                seq = list(seq)
                try:
                    seq[key] = val


                except (TypeError):
                    pass
    seq=''.join(seq)
    g.add_line("S\tNEW\t*")
    f = g.segments
    for i in f:
        if i == "NEW":
            print("check LEN")
            new_line = i
            new_line.name = str(edge) + "_" + str(clN)
            new_line.sid = str(edge) + "_" + str(clN)
            new_line.sequence = seq
            new_line.dp = 35  # coverage
            print("edge added:"+str(new_line.name))
            print(len(new_line.sequence))

def change_link(old, new, neighbor, neighbor_cl=None):  #смениить на addlink
    if old==neighbor:
        exit()
    print("CNANGE I!!! "+old+" "+neighbor)
    old=str(old)+"\t"
    new=str(new)+"\t"

    for i in g.edges:
        if re.search(old, str(i)):

            if neighbor_cl == None or g.segment("%s_%s\t" % (neighbor, neighbor_cl))==None:

                if re.search(neighbor, str(i)):
                    print("1")
                    print(str(i).replace(old, new))
                    print(i)

                    g.add_line(str(i).replace(old, new))
                    #g.rm(i)
            else:
                if re.search(neighbor, str(i)):
                    print("2")
                    print(neighbor)


                    print(str(i).replace(old, new).replace(str(neighbor) + "\t", "%s_%s\t" % (neighbor, neighbor_cl)))
                    print(i)

                    g.add_line(str(i).replace(old, new).replace(str(neighbor) + "\t", "%s_%s\t" % (neighbor, neighbor_cl)))
                    #g.rm(i)

def remove_link(edge, neighbor):
    print("remove "+edge+" "+ neighbor)
    edge=str(edge)+"\t"
    for i in g.edges:
        if re.search(edge, str(i)):
            if re.search(neighbor, str(i)):
                print(i)
                g.rm(i)



import gfapy
g = gfapy.Gfa.from_file("/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph.gfa")
I = 1000
AF = 0.1
full_cl={}

print("CREATING NEW EDGES")
for edge in edges:
    print(edge)
    #cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
        snp=None
        SNP_pos=read_snp(snp)
        data=read_bam(bam,SNP_pos)
        full_clusters = []
        ln = int(pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4])
        clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
        for cluster in clusters:

            clStart = int(ln)
            clStop = 0
            for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
                start = int(data[read]["Start"])
                stop = int(data[read]["Stop"])
                if start < clStart:
                    clStart = start
                if stop > clStop:
                    clStop = stop

            if clStart == 0 and clStop == ln - 1:
                full_clusters.append(cluster)
                add_child_edge(edge, cluster, g, cl, SNP_pos, data)

        full_cl[edge] = full_clusters

    except(FileNotFoundError): continue



print(full_cl)
print(len(full_cl))




print("CREATING NEW LINKS")
for edge in full_cl.keys():

    clusters=full_cl[edge]
    try:
        cl=pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
    except(FileNotFoundError): pass
    gaf=pd.read_csv("~/MT/gaf3noalt/gaf_%s.csv"  % (edge))
    new_links=[]
    #print(gaf)
    for clN in clusters:
        print("%s_%s" % (edge, clN))
        reads=list(cl.loc[cl['Cluster'] == str(clN), 'ReadName'])
        als =gaf.loc[gaf['ReadName'].isin(reads)]
        als=als.to_dict('split')['data']
        head={}
        tail={}


        for aln in als:
            #print(aln)
            al=aln[2]
            read=aln[1]
            if al.split("%s" % (edge))[0]=='>':
                try:tail[read]=al.split(">")[2]
                except(IndexError): pass
            elif al.split("%s" % (edge))[0]=='<':
                try:
                    head[read]=al.split("<")[2]
                except(IndexError): pass

            elif len(al.split("%s" % (edge)))>1:
                try: tail[read]=al.split("<")[1]
                except(IndexError): pass
            elif len(al.split(">%s" % (edge)))>1:
                try: head[read]=al.split(">")[len(al.split(">"))-2]
                except(IndexError): pass
        #print(set(head.values()))
        #print(tail)


    #ищем связи начала
    #в каком кластере эти риды?
        for h in set(head.values()):

            try: cl_h = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (h, I, AF), keep_default_na=False)
            except(FileNotFoundError):
                break
            reads=[]
            for k, v in head.items():
                if v==h:
                    reads.append(k)
            h_cl = cl_h.loc[cl_h['ReadName'].isin(reads), 'Cluster']
            h_cl=list(set(h_cl))
            #print(h_cl)
        #для каждого h_cl
        #проверить есть ли такое ребро
        #print(h_cl)
            if len(h_cl)==0:
                print("no child")
                try: #g.add_line("L\t%s\t+\t%s_%s\t-\t*" % (h, edge, clN))
                    change_link(edge,"%s_%s" % (edge, clN), h)
                    new_links.append(h)


                except:
                    pass
                print("addd links to " + str("%s" % (h)))
            for i in h_cl:
                if '%s_%s' %(h, i) not in g.segment_names:
                    print("no child")
                    try: #g.add_line("L\t%s\t+\t%s_%s\t-\t*" % (h, edge, clN))
                        change_link(edge, "%s_%s" % (edge, clN), h)
                        new_links.append(h)
                    except:pass
                    print("addd links to " + str("%s" % (h)))
                #нет подребра, добавляем связь к головной ноде
                    pass
                else:
            #ПРОВЕРИТЬ НАПРАВЛЕНИЕ
            #есть подребро==>cделать связь edge_clN c h, h_cl[0]
                    try:  #g.add_line("L\t%s_%s\t+\t%s_%s\t-\t*"  %(h, i,edge, clN))
                        change_link(edge, "%s_%s" % (edge, clN), h, i)
                        new_links.append(h)
                    except:pass
                    print("addd links to "+str("%s_%s" %(h, i)))

                #удалить связь edge_clN c h
    # ищем связи конца
        for t in set(tail.values()):
            try:cl_t = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (h, I, AF), keep_default_na=False)
            except(FileNotFoundError): break
            reads = []
            for k, v in head.items():
                if v == h:
                    reads.append(k)
            t_cl = cl_t.loc[cl_t['ReadName'].isin(reads), 'Cluster']
            t_cl = list(set(h_cl))
        # print(h_cl)
        # для каждого h_cl
        # проверить есть ли такое ребро
        # print(h_cl)
            if len(t_cl) == 0:
                print("no child")
                try: #g.add_line("L\t%s\t-\t%s_%s\t+\t*" % (t, edge, clN))
                    change_link(edge, "%s_%s" % (edge, clN), t)
                    new_links.append(t)
                except:
                    pass
                print("addd links to " + str("%s" % (h)))
            for i in t_cl:
                if '%s_%s' % (h, i) not in g.segment_names:
                    print("no child")
                    try: #g.add_line("L\t%s\t-\t%s_%s\t+\t*" % (t, edge, clN))
                        change_link(edge, "%s_%s" % (edge, clN), t)
                        new_links.append(t)
                    except:
                        pass
                    print("addd links to " + str("%s" % (t)))
                # нет подребра, добавляем связь к головной ноде
                    pass
                else:
                # ПРОВЕРИТЬ НАПРАВЛЕНИЕ
                # есть подребро==>cделать связь edge_clN c h, h_cl[0]
                    try:#g.add_line("L\t%s_%s\t-\t%s_%s\t+\t*" % (t, i, edge, clN))
                        change_link(edge, "%s_%s" % (edge, clN), t,i)
                        new_links.append(t)
                    except:pass
                    print("addd links to " + str("%s_%s" % (t, i)))
    for i in new_links:
        remove_link(edge,i)

    #g.rm(edge)  #временно удаляю edge

    gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph_transformed.gfa")

#g.add_line("L\tedge_%s_%s\t+\tedge_263\t-\t*" %(149, 1))
#line = g.segment("L\tedge_1\t+\tedge_267\t-\t*")
#print(line)



#print(g.edges)
#g.rm("L\tedge_3\t+\tedge_7\t+\t*")
#поменять cov старой

#g.RECORDS_WITH_NAME





#print(g.segment_names)
#print(g.try_get_segment('edge_100'))
print(g.segment('edge_100_1'))




gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/flye_3ecoli_sim/assembly_graph_transformed.gfa")
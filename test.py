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
from cluster_postprocess import build_adj_matrix_clusters
import sys, os, subprocess
import pygraphviz as gv
import re
import gfapy
from collections import Counter
from collections import deque

# import pylab

edges = ['edge_1', 'edge_10', 'edge_100', 'edge_101', 'edge_102', 'edge_103', 'edge_104', 'edge_105', 'edge_106',
         'edge_107', 'edge_108', 'edge_109', 'edge_11', 'edge_110', 'edge_111', 'edge_112', 'edge_113', 'edge_114',
         'edge_115', 'edge_116', 'edge_12', 'edge_13', 'edge_14', 'edge_15', 'edge_16', 'edge_17', 'edge_18', 'edge_19',
         'edge_2', 'edge_20', 'edge_21', 'edge_22', 'edge_23', 'edge_24', 'edge_25', 'edge_26', 'edge_27', 'edge_28',
         'edge_29', 'edge_3', 'edge_30', 'edge_31', 'edge_32', 'edge_33', 'edge_34', 'edge_35', 'edge_36', 'edge_37',
         'edge_38', 'edge_39', 'edge_4', 'edge_40', 'edge_41', 'edge_42', 'edge_43', 'edge_44', 'edge_45', 'edge_46',
         'edge_47', 'edge_48', 'edge_49', 'edge_5', 'edge_50', 'edge_51', 'edge_52', 'edge_53', 'edge_54', 'edge_55',
         'edge_56', 'edge_57', 'edge_58', 'edge_59', 'edge_6', 'edge_60', 'edge_61', 'edge_62', 'edge_63', 'edge_64',
         'edge_65', 'edge_66', 'edge_67', 'edge_68', 'edge_69', 'edge_7', 'edge_70', 'edge_71', 'edge_72', 'edge_73',
         'edge_74', 'edge_75', 'edge_76', 'edge_77', 'edge_78', 'edge_79', 'edge_8', 'edge_80', 'edge_81', 'edge_82',
         'edge_83', 'edge_84', 'edge_85', 'edge_86', 'edge_87', 'edge_88', 'edge_89', 'edge_9', 'edge_90', 'edge_91',
         'edge_92', 'edge_93', 'edge_94', 'edge_95', 'edge_96', 'edge_97', 'edge_98', 'edge_99']
# edges=[ 'edge_10', 'edge_101',  'edge_103', 'edge_104', 'edge_105', 'edge_108', 'edge_11','edge_110', 'edge_111',  'edge_114', 'edge_12', 'edge_15',  'edge_2', 'edge_21', 'edge_22', 'edge_23','edge_28', 'edge_29', 'edge_31',  'edge_35', 'edge_36','edge_4','edge_42', 'edge_47', 'edge_48',  'edge_5','edge_54','edge_60',  'edge_62',  'edge_65', 'edge_71', 'edge_72',  'edge_74', 'edge_75', 'edge_79', 'edge_84',  'edge_88', 'edge_89', 'edge_90', 'edge_91', 'edge_92', 'edge_94', 'edge_96', 'edge_97', 'edge_98', 'edge_99']


# edges=['edge_96', "edge_43"]

# edges = ['edge_99', 'edge_96', 'edge_98']


edges = ['edge_20']

R = 1
I = 1000
clipp = 100
min_mapping_quality = 20
min_base_quality = 0
min_al_len = 1000
de_max = 0.05
AF = 0.1

bam = "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/sim.3ecoliv2.bam"
stats = open('output/stats_clusters.txt', 'a')
stats.write("Edge" + "\t" + "Fill Clusters" + "\t" + "Full Paths Clusters" + "\n")
stats.close()


def read_snp(snp, edge):
    SNP_pos = []
    if snp == None:
        snpos = 'bcftools mpileup -r {} {} --no-reference -I --no-version --annotate FORMAT/AD | bcftools query -f  "%CHROM %POS [ %AD %DP]\n" >output/vcf/vcf_{}.txt'.format(
            edge, bam, edge)
        subprocess.check_output(snpos, shell=True, capture_output=False)
        with open("output/vcf/vcf_%s.txt" % edge) as f:
            lines = f.readlines()
            for line in lines:
                try:
                    AlFreq = int(str(line.split()[2]).split(',')[2]) / int(line.split()[3])
                except(IndexError):
                    AlFreq = 0
                if AlFreq > AF:
                    SNP_pos.append(line.split()[1])

    else:
        vcf = open(snp, "rt")
        for line in vcf:
            if line.split()[0] == edge:
                SNP_pos.append(line.split()[1])
    # print(str(len(SNP_pos)) + " SNPs found")
    # print(SNP_pos)
    return (SNP_pos)


def read_bam(bam, edge, SNP_pos, clipp=clipp):
    print("read bam")
    print(edge)
    bamfile = pysam.AlignmentFile(bam, "rb")
    data = {}
    ln = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4]
    for pos in SNP_pos:
        for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools', min_base_quality=0,
                                           ignore_overlaps=False, min_mapping_quality=20,
                                           ignore_orphans=False, truncate=True):

            for pileupread in pileupcolumn.pileups:
                clipping = 0
                start = pileupread.alignment.get_reference_positions()[0]
                stop = pileupread.alignment.get_reference_positions()[
                    len(pileupread.alignment.get_reference_positions()) - 1]

                de = float(str(pileupread).split('\t')[11].split('), (')[8].split(',')[1])
                for i in pileupread.alignment.cigartuples:
                    if i[0] == 4 or i[0] == 5:
                        if i[1] > clipp:
                            clipping = 1

                if (clipping == 1 or (stop - start) < min_al_len or de > de_max) and (
                        int(start) != 0 and int(stop) != int(ln) - 1):
                    continue


                else:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        try:
                            data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[
                                pileupread.query_position]

                        except (KeyError):
                            data[pileupread.alignment.query_name] = {}
                            data[pileupread.alignment.query_name]["Start"] = \
                                pileupread.alignment.get_reference_positions()[
                                    0]
                            data[pileupread.alignment.query_name]["Stop"] = \
                                pileupread.alignment.get_reference_positions()[
                                    len(pileupread.alignment.get_reference_positions()) - 1]

                            data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[
                                pileupread.query_position]

    bamfile.close()
    return (data)


def add_child_edge(edge, clN, g, cl, SNP_pos, data, left, righ):
    f = g.segments
    seq=[]
    # dat = {0: "A", 1: "A", 4: "A"}
    cons = build_data(cl, SNP_pos, data)
    cl_consensuns = cluster_consensuns(cl, clN, SNP_pos, data, cons)
    print(cl_consensuns[clN])
    for i in f:
        if i == edge:
            seq = i.sequence
            seq = list(seq)
            for key, val in cl_consensuns[clN].items():
                try:
                    seq[int(key)] = val

                except (ValueError):
                    continue

    # print(len(seq))


    seq = ''.join(seq)

    seq = seq[left:righ]

    g.add_line("S\tNEW\t*")
    f = g.segments
    for i in f:
        if i == "NEW":
            new_line = i
            new_line.name = str(edge) + "_" + str(clN)
            new_line.sid = str(edge) + "_" + str(clN)
            new_line.sequence = seq
            new_line.dp = 35  # coverage
            print("edge added:" + str(new_line.name))
            print(left, righ)
    # gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph_transformed.gfa")


def remove_link(edge, neighbor):
    # print("remove "+edge+" "+ neighbor)
    edge = str(edge) + "\t"
    for i in g.edges:
        if re.search(edge, str(i)):
            if re.search(neighbor, str(i)):
                print(i)
                g.rm(i)


def build_paths_graph(edge, data, SNP_pos, cl):
    cons = build_data(cl, SNP_pos, data)
    M = build_adj_matrix_clusters(cons, SNP_pos, cl)
    M = change_w(M, 1)
    G_vis = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    for node in G_vis.nodes():
        neighbors = nx.all_neighbors(G_vis, node)
        for neighbor in list(neighbors):
            n_paths = nx.algorithms.all_simple_paths(G_vis, node, neighbor)
            for n_path in list(n_paths):
                if len(n_path) == 3:
                    print(n_path)
                    G_vis.remove_edge(n_path[0], n_path[1])

    G_vis = nx.nx_agraph.to_agraph(G_vis)
    G_vis.layout(prog="neato")
    G_vis.draw("output/graphs/full_paths_cluster_GV_graph_%s.png" % (edge))
    G = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
    return (G)


def find_full_paths(G, paths_roots, paths_leafs):
    for node in G.nodes():
        neighbors = nx.all_neighbors(G, node)
        for neighbor in list(neighbors):
            n_paths = nx.algorithms.all_simple_paths(G, node, neighbor)
            for n_path in list(n_paths):
                if len(n_path) == 3:
                    print(n_path)
                    G.remove_edge(n_path[0], n_path[1])
    paths = []
    for root in paths_roots:
        paths_nx = nx.algorithms.all_simple_paths(G, root, paths_leafs)
        for path in list(paths_nx):
            paths.append(path)
    return (paths)


def add_link(first, second):
    str = 'L	first_edge	+	second_edge	+	0M	RC:i:42'
    try:
        g.add_line(str.replace('first_edge', first).replace('second_edge', second))
        print("line added " + first + " " + second)
    except(gfapy.NotUniqueError): pass
        #print("dd")


def add_path_links(edge, paths):
    str = 'L	first_edge	+	second_edge	+	0M	RC:i:42'
    for path in paths:

        for i in range(0, len(path) - 1):
            try:
                print("path added " + str(path[i]) + " " + str(path[i + 1]))
            except:
                print("path added")
            try:
                g.add_line(str.replace('first_edge', "%s_%s" % (edge, path[i])).replace('second_edge',
                                                                                        "%s_%s" % (edge, path[i + 1])))

            except(gfapy.NotUniqueError):
                continue


def add_path_edges ( edge,g,cl, data, SNP_pos, ln, paths, G,paths_roots,paths_leafs, clusters_borders):
    path_cl = []
    for path in paths[edge]:
        for member in path:
            path_cl.append(member)
    cut_l = {}
    cut_r = {}
    for path_cluster in set(path_cl):
        if path_cluster in paths_roots:
            cut_l[path_cluster] = 0
            cut_r[path_cluster] = None
        elif path_cluster in paths_leafs:
            cut_r[path_cluster] = ln - 1
            cut_l[path_cluster] = None
        else:
            cut_l[path_cluster] = None
            cut_r[path_cluster] = None

    while None in cut_l.values():
        for member in cut_l.keys():
            if cut_l[member] != None and cut_r[member] == None:
                # print("NEW")
                # print(member)
                Q = deque()
                L = []
                R = []
                for path in paths[edge]:
                    try:
                        L.append(path[path.index(member) + 1])
                        Q.append(path[path.index(member) + 1])
                    except (ValueError, IndexError):
                        continue

                visited = []
                while Q:
                    n = Q.pop()
                    visited.append(n)
                    if n in L:
                        for path in paths[edge]:
                            try:
                                R.append(path[path.index(n) - 1])
                                if path[path.index(n) - 1] not in visited:
                                    Q.append(path[path.index(n) - 1])


                            except (ValueError, IndexError):
                                continue
                    else:
                        for path in paths[edge]:
                            try:
                                L.append(path[path.index(n) + 1])
                                if path[path.index(n) + 1] not in visited:
                                    Q.append(path[path.index(n) + 1])

                            except (ValueError, IndexError):
                                continue
                l_borders = []
                r_borders = []
                for i in L:
                    l_borders.append(int(clusters_borders[i][0]))
                for i in R:
                    r_borders.append(int(clusters_borders[i][1]))

                border = max(l_borders) + (min(r_borders) - max(l_borders)) // 2
                for i in L:
                    cut_l[i] = border

                for i in R:
                    cut_r[i] = border


    for path_cluster in set(path_cl):
        add_child_edge(edge, path_cluster, g,  cl, SNP_pos, data, cut_l[path_cluster], cut_r[path_cluster])
        print("edge added ask")
        G.remove_node(path_cluster)
    return(path_cl)


def add_parent_subedge(edge, g, left, right):
    f = g.segments
    for i in f:
        if i == edge:
            seq = i.sequence
            seq = list(seq)
    seq = ''.join(seq)
    seq = seq[left:right]

    g.add_line("S\tNEW\t*")
    f = g.segments
    for i in f:
        if i == "NEW":
            new_line = i
            new_line.name = "%s[%s,%s]" % (edge, left, right)
            new_line.sid = "%s[%s,%s]" % (edge, left, right)
            new_line.sequence = seq
            new_line.dp = 35  # coverage
            print("edge added:" + str(new_line.name))
            print(left, right)


def add_parent_subunits(edge, g, parent_borders, ln):
    parent_subunits = []
    parent_borders.append(0)
    parent_borders.append(ln - 1)
    parent_borders=list(set(parent_borders))
    subunits_borderline[edge]= [[], []]



    for i in range(0, len(parent_borders)-1):
        add_parent_subedge(edge, g, sorted(parent_borders)[i], sorted(parent_borders)[i + 1])
        parent_subunits.append("%s[%s,%s]" % (edge, sorted(parent_borders)[i], sorted(parent_borders)[i + 1]))

        if parent_borders[i]==0:
            subunits_borderline[edge][0].append("%s[%s,%s]" % (edge, sorted(parent_borders)[i], sorted(parent_borders)[i + 1]))
        if parent_borders[i+1]==ln-1:
            subunits_borderline[edge][1].append("%s[%s,%s]" % (edge, sorted(parent_borders)[i], sorted(parent_borders)[i + 1]))
    return (parent_subunits)


g = gfapy.Gfa.from_file(
    "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph.gfa")
I = 1000
AF = 0.1
full_cl = {}
full_paths = {}
paths = {}
full_path_clusters = {}
subunits_borderline={}
connected_subunits={}


def transform_graph(i):
    edge = edges[i]
    print(edge)
    full_paths_roots = []
    full_paths_leafs = []
    full_clusters = []
    paths_roots = []
    paths_leafs = []
    part_clusters = []
    edges_to_remove=[]
    # CHECK if FULL
    # cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
        snp = None
        SNP_pos = read_snp(snp, edge)
        data = read_bam(bam, edge, SNP_pos)
        ln = int(pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4])
        clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))

        clusters_borders = {}
        G = build_paths_graph(edge, data, SNP_pos, cl)

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
                add_child_edge(edge, cluster, g, cl, SNP_pos, data, 0, ln - 1)
                G.remove_node(cluster)

            elif clStart == 0:
                full_paths_roots.append(cluster)

            elif clStop == ln - 1:
                full_paths_leafs.append(cluster)

            clusters_borders[cluster] = [clStart, clStop]
        full_cl[edge] = full_clusters
        # CHECK if SEMI-FULL

        try:

            full_paths[edge] = find_full_paths(G,full_paths_roots, full_paths_leafs)
        except(ValueError):
            pass
        # continue
        path_cl = add_path_edges(edge,g,cl, data, SNP_pos, ln,full_paths, G,full_paths_roots, full_paths_leafs,clusters_borders)
        add_path_links(edge, full_paths[edge])

        gfapy.Gfa.to_file(g,
                          "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph_transformed.gfa")
        full_path_clusters[edge] = set(path_cl)
        print(len(set(path_cl)))
        '''
        lis = list(nx.topological_sort(nx.line_graph(G)))

        first = []
        last = []
        parent_borders = []
        for i in lis:
            first.append(i[0])
            last.append(i[1])



        for i in lis:

            if last.count(i[0]) == 0:
                paths_roots.append(i[0])
                parent_borders.append(clusters_borders[i[0]][0])
            if first.count(i[1]) == 0:
                paths_leafs.append(i[1])
                parent_borders.append(clusters_borders[i[1]][1])
        try:
            paths[edge] = find_full_paths(G, paths_roots, paths_leafs)

        except(ValueError):
            pass


        print(paths[edge])

        subclusters=list(nx.isolates(G))
        for i in subclusters:
            parent_borders.append(clusters_borders[i][0])
            parent_borders.append(clusters_borders[i][1])
        print(parent_borders)
        for subcluster in set(subclusters):
            add_child_edge(edge, subcluster, g, cl, SNP_pos, data, clusters_borders[subcluster][0], clusters_borders[subcluster][1])
        print(paths)
        path_cl=add_path_edges(edge,g,cl, data, SNP_pos, ln,paths, G, paths_roots, paths_leafs,clusters_borders)

        add_path_links(edge, paths[edge])

        parent_subunits=add_parent_subunits(edge, g, parent_borders, ln)
        if len (parent_borders)>0:
            for i in paths_roots:
                for subedge in parent_subunits:
                    if re.search("%s.*,%s\]" % (edge, clusters_borders[i][0]), subedge):
                        add_link(subedge, str(i))
            for i in paths_leafs:
                for subedge in parent_subunits:
                    if re.search("%s\[%s,.*]" % (edge, clusters_borders[i][1]), subedge):
                        add_link(str(i), subedge)

            connected_subunits[edge]=[]
            for i in subclusters:
                for subedge in parent_subunits:
                    if re.search("%s\[%s,.*]" % (edge, clusters_borders[i][1]), subedge):
                        add_link('%s_%s' % (edge,str(i)), subedge)
                        connected_subunits[edge].append(subedge)
                    if re.search("%s.*,%s\]" % (edge, clusters_borders[i][0]), subedge):
                        add_link(subedge, '%s_%s' % (edge,str(i)))
                        connected_subunits[edge].append(subedge)
            print("connect SUBUNITS")
            print(connected_subunits[edge])
            for connected_subunit1 in connected_subunits[edge]:

                for connected_subunit2 in connected_subunits[edge]:
                    if re.search("%s\[%s,.*]" % (edge, re.sub('\]',"",re.sub("%s" % edge, "", connected_subunit1).split(',')[1])), connected_subunit2):
                        add_link(connected_subunit1, connected_subunit2)
                    if re.search("%s.*,%s\]" % (edge, re.sub('\[',"",re.sub("%s" % edge, "", connected_subunit1).split(',')[0])), connected_subunit2):
                        add_link(connected_subunit2, connected_subunit1)


            for i in parent_subunits:
                for seg in g.segments:
                    if seg.name==i and i not in connected_subunits[edge]:
                        g.rm(seg)
                        #print("REMOVE")
    '''
    except(FileNotFoundError):
        pass
    # continue


    stats = open('output/stats_clusters.txt', 'a')
    fcN = 0
    fpN = 0
    try:
        fcN = len(full_cl[edge])
    except(KeyError):
        pass
    try:
        fpN = len(full_paths[edge])
    except(KeyError):
        pass

    othercl=len(clusters)-fcN-fpN

    if othercl==0:
        edges_to_remove.append(edge)
    stats.write(edge + "\t" + str(fcN) + "\t" + str(fpN) + "\t" + str(othercl) +"\n")
    stats.close()

    # print(full_cl)
    # print(len(full_cl))
    # print(full_paths)
    # print(full_paths_roots,full_paths_leafs)

    print("CREATING NEW LINKS")
    # for edge in edges:
    clusters = []
    print()
    try:
        clusters = full_cl[edge]
    except(KeyError):pass

    try:
        clusters = clusters + list(full_path_clusters[edge])
    except:
        pass

    try:
        clusters = clusters +list(paths_leafs)+list(paths_roots)
    except:
        pass
    # continue
    print(clusters)

    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
    except(FileNotFoundError):
        pass
    # gaf = pd.read_csv("~/MT/gaf3noalt/gaf_%s.csv" % (edge))
    gaf = pd.read_csv("~/MT/gaf3noalt/gaf_common.csv")
    new_links = []
    for clN in clusters:
        print("%s_%s" % (edge, clN))
        reads = list(cl.loc[cl['Cluster'] == clN, 'ReadName'])
        als = gaf.loc[gaf['ReadName'].isin(reads)]
        als = als.to_dict('split')['data']
        head = {}
        tail = {}
        for aln in als:
            al = aln[2]
            read = aln[1]
            if re.search(".*%s<" % edge, al):
                h = re.sub("<.*", "",
                           re.sub(">.*", "", re.sub(".*%s<" % edge, "", al, count=0, flags=0), count=0, flags=0),
                           count=0,
                           flags=0)
                if len(h) > 0:
                    head[read] = h

            if re.search(".*>%s.*" % edge, al):
                h = re.sub(".*>", "",
                           re.sub(".*<", "", re.sub(">%s.*" % edge, "", al, count=0, flags=0), count=0, flags=0),
                           count=0,
                           flags=0)
                if len(h) > 0:
                    head[read] = h

            if re.search(".*%s>" % edge, al):
                t = re.sub("<.*", "",
                           re.sub(">.*", "", re.sub(".*%s>" % edge, "", al, count=0, flags=0), count=0, flags=0),
                           count=0,
                           flags=0)
                if len(t) > 0:
                    tail[read] = t
            if re.search(".*<%s.*" % edge, al):
                t = re.sub(".*>", "",
                           re.sub(".*<", "", re.sub("<%s.*" % edge, "", al, count=0, flags=0), count=0, flags=0),
                           count=0,
                           flags=0)
                if len(t) > 0:
                    tail[read] = t

        for h in set({k for k, v in Counter(head.values()).items() if v > 2}):
            try:
                cl_h = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (h, I, AF), keep_default_na=False)
            except(FileNotFoundError):
                # print("break")
                add_link(h, "%s_%s" % (edge, clN))
                new_links.append(h)
                break
            reads = []
            for k, v in head.items():
                if v == h:
                    reads.append(k)
            h_cl = cl_h.loc[cl_h['ReadName'].isin(reads), 'Cluster']
            h_cl = list(set(h_cl))
            print(h_cl)

            if len(h_cl) == 0:
                # print("no child")
                try:  # g.add_line("L\t%s\t+\t%s_%s\t-\t*" % (h, edge, clN))
                    # change_link(edge,"%s_%s" % (edge, clN), h)
                    add_link(h, "%s_%s" % (edge, clN))
                    new_links.append(h)


                except:
                    pass
                print("addd links to " + str("%s" % (h)))
            for i in h_cl:
                if '%s_%s' % (h, i) not in g.segment_names:
                    # print("no child")
                    try:  # g.add_line("L\t%s\t+\t%s_%s\t-\t*" % (h, edge, clN))
                        add_link(h, "%s_%s" % (edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), h)
                        new_links.append(h)
                    except:
                        pass
                    print("addd links to " + str("%s" % (h)))
                    # нет подребра, добавляем связь к головной ноде
                    pass
                else:
                    try:  # g.add_line("L\t%s_%s\t+\t%s_%s\t-\t*"  %(h, i,edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), h, i)
                        add_link("%s_%s" % (h, i), "%s_%s" % (edge, clN))
                        new_links.append(h)
                    except:
                        pass
                    print("addd links to " + str("%s_%s" % (h, i)))

        for t in set({k for k, v in Counter(tail.values()).items() if v > 2}):
            try:
                cl_t = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (t, I, AF), keep_default_na=False)
            except(FileNotFoundError):
                add_link("%s_%s" % (edge, clN), t)
                new_links.append(t)
                break
            reads = []
            for k, v in head.items():
                if v == h:
                    reads.append(k)
            t_cl = cl_t.loc[cl_t['ReadName'].isin(reads), 'Cluster']
            t_cl = list(set(t_cl))

            if len(t_cl) == 0:
                print("no child")
                try:  # g.add_line("L\t%s\t-\t%s_%s\t+\t*" % (t, edge, clN))
                    # change_link(edge, "%s_%s" % (edge, clN), t)
                    add_link("%s_%s" % (edge, clN), t)
                    new_links.append(t)
                except:
                    pass
                print("addd links to " + str("%s" % (t)))
            for i in t_cl:
                if '%s_%s' % (h, i) not in g.segment_names:
                    # print("no child")
                    try:  # g.add_line("L\t%s\t-\t%s_%s\t+\t*" % (t, edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), t)
                        add_link("%s_%s" % (edge, clN), t)
                        new_links.append(t)
                    except:
                        pass
                    print("addd links to " + str("%s" % (t)))

                    pass
                else:

                    try:  # g.add_line("L\t%s_%s\t-\t%s_%s\t+\t*" % (t, i, edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), t,i)
                        add_link("%s_%s" % (edge, clN), "%s_%s" % (t, i))
                        new_links.append(t)
                    except:
                        pass
                    print("addd links to " + str("%s_%s" % (t, i)))
        for ed in g.edges:
           if ed.name in edges_to_remove:
               g.rm(ed)
        #for i in new_links:
            #remove_link(edge, i)
        gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph_transformed.gfa")


for i in range(0, len(edges)):
   transform_graph(i)

'''
for edge in edges:
   ###здесь еще нужно проверять направление ребра!!!= или -
    print(subunits_borderline)
    for ed in g.edges:
        if re.search("%s\t" %edge, str(ed)):
            if str(ed).split("\t")[1]==edge:
                for i in subunits_borderline[edge][1]:
                    if i in connected_subunits[edge]:
                        g.add_line(str(ed).replace("%s\t" %edge, "%s\t" % i))
                        g.rm(ed)
                        print(str(ed).replace("%s\t" %edge, "%s\t" % i))
            elif str(ed).split("\t")[3]==edge:
                for i in subunits_borderline[edge][0]:
                    if i in connected_subunits[edge]:
                        g.add_line(str(ed).replace("%s\t" %edge, "%s\t" % i))
                        g.rm(ed)
                        print(str(ed).replace("%s\t" % edge, "%s\t" % i))

        gfapy.Gfa.to_file(g,
                          "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph_transformed.gfa")

'''


import multiprocessing


def parallel(edges):
    pool = multiprocessing.Pool(1)
    pool.map(transform_graph, range(0, len(edges)))
    pool.close()


#if __name__ == "__main__":
    #parallel(edges)

# REMOVE LINKS


gfapy.Gfa.to_file(g,
                  "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph_transformed.gfa")

# g.add_line("L\tedge_%s_%s\t+\tedge_263\t-\t*" %(149, 1))
# line = g.segment("L\tedge_1\t+\tedge_267\t-\t*")
# print(line)


# print(g.edges)
# g.rm("L\tedge_3\t+\tedge_7\t+\t*")
# поменять cov старой

# g.RECORDS_WITH_NAME


# print(g.segment_names)
# print(g.try_get_segment('edge_100'))
# print(g.segment('edge_100_1'))

'''
def change_link(old, new, neighbor, neighbor_cl=None):  # смениить на addlink
    if old == neighbor:
        exit()
    # print("CNANGE I!!! "+old+" "+neighbor)
    old = str(old) + "\t"
    new = str(new) + "\t"
    temp = []

        for i in g.edges:
        if re.search(old, str(i)):

            if neighbor_cl == None or g.segment("%s_%s\t" % (neighbor, neighbor_cl)) == None:

                if re.search(neighbor, str(i)):
                    print("1")
                    print(str(i).replace(old, new))
                    print(i)

                    g.add_line(str(i).replace(old, new))
                    # g.rm(i)
            else:
                if re.search(neighbor, str(i)):
                    print("2")
                    print(neighbor)

                    print(str(i).replace(old, new).replace(str(neighbor) + "\t", "%s_%s\t" % (neighbor, neighbor_cl)))
                    print(i)

                    g.add_line(
                        str(i).replace(old, new).replace(str(neighbor) + "\t", "%s_%s\t" % (neighbor, neighbor_cl)))
                    # g.rm(i)'''

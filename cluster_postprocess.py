import csv
import pysam
import pandas as pd
from collections import Counter
import sys,os,subprocess
import networkx as nx
import matplotlib.pyplot as plt
from build_adj_matrix import *
import pygraphviz as gv
import pylab
from community_detection import find_communities


def build_data_cons(cl,SNP_pos, data):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
    #clusters = sorted(set(cl['Cluster'].values))
    cons = {}
    for cluster in clusters:
        cons=cluster_consensuns(cl,cluster,SNP_pos, data, cons)
    #print(cons)
    return(cons)


def cluster_consensuns(cl,cluster,SNP_pos, data, cons):
    strange = 0
    val = {}
    clSNP = []
    for pos in SNP_pos:
        # print(pos)
        npos = []
        for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
            # print(read)
            try:
                npos.append(data[read][pos])
            except(KeyError):
                continue
        # print(npos)

        try:
            if len(npos) >= 2:
                # print(Counter(npos).most_common())[1][1]
                if int(Counter(npos).most_common()[0][1]) >= 2:
                    val[pos] = Counter(npos).most_common()[0][0]
            if int(Counter(npos).most_common()[1][1]) >= 2:
                strange = 1
                clSNP.append(pos)
        except(IndexError):
            continue

    if strange == 1:
        val["Strange"] = 1
    else:
        val["Strange"] = 0
    val["clSNP"] = clSNP

    clStart = 1000000000000  # change fo ln
    clStop = 0
    clCov=0

    for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
        start = int(data[read]["Start"])
        stop = int(data[read]["Stop"])
        clCov=clCov+(stop-start)
        if start < clStart:
            clStart = start
        if stop > clStop:
            clStop = stop

    clCov=clCov/(clStop-clStart)
    val["Stop"] = clStop
    val["Start"] = clStart
    val["Cov"]=clCov
    cons[cluster] = val
    return (cons)


def split_cluster(cl,cluster, data,clSNP, bam, edge, child_clusters, R, I):
    print("Strange cluster detected")
    reads=sorted(set(cl.loc[cl['Cluster'] == cluster]['ReadName'].values))
    if cluster==1000000:
        print("Build na matrix")
        m = build_adj_matrix2(cl[cl['Cluster'] == cluster], data, clSNP, I, bam, edge, R, only_with_common_snip=False)
    else:
        m=build_adj_matrix2(cl[cl['Cluster'] == cluster], data, clSNP, I, bam,edge,R)
    m = remove_edges(m, 1)
    m.columns=range(0,len(cl[cl['Cluster'] == cluster]['ReadName']))
    m.index=range(0,len(cl[cl['Cluster'] == cluster]['ReadName']))
    m = change_w(m, R)
    G = nx.from_pandas_adjacency(m)
    #nx.draw(G, nodelist=G.nodes(), with_labels=False, width=0.03, node_size=3, font_size=5)
    #plt.savefig("output/test.png" , format="PNG", dpi=300)
    #plt.close()
    cluster_membership = find_communities(G)
    #print(cluster_membership)
    clN=0
    uncl=0
    reads = cl[cl['Cluster'] == cluster]['ReadName'].values

    for value in set(cluster_membership.values()):
        group = [k for k, v in cluster_membership.items() if v == value]

        if len(group) > 2:
            clN = clN + 1
            child_clusters.append(cluster+10000+clN)
            for i in group:
                cl.loc[cl['ReadName'] == reads[i], "Cluster"] =  cluster+10000+clN

        else:
            uncl = uncl + 1
            for i in group:
                if cluster==1000000:
                    cl.loc[cl['ReadName'] == reads[i], "Cluster"] = 'NA'
                else:
                    cl.loc[cl['ReadName'] == reads[i], "Cluster"] = 1000000

    print(str(clN)+" new clusters found")


def distance_clusters(first_cl,second_cl, cons,SNP_pos, only_with_common_snip=True):
    d=-1
    #print("distance "+str(first_cl)+" "+str(second_cl))
    firstSNPs=list(cons[first_cl].keys())
    firstSNPs.remove('clSNP')
    firstSNPs.remove('Strange')
    firstSNPs.remove('Stop')
    firstSNPs.remove('Start')
    firstSNPs.remove('Cov')
    secondSNPs=list(cons[second_cl].keys())
    secondSNPs.remove('clSNP')
    secondSNPs.remove('Strange')
    secondSNPs.remove('Stop')
    secondSNPs.remove('Start')
    secondSNPs.remove('Cov')
    commonSNP=sorted(set(firstSNPs).intersection(secondSNPs))
    #print(commonSNP)
    try:
        #if abs(int(commonSNP[len(commonSNP)-1])-int(commonSNP[0]))>=500:
        #if 1000 >= 500:
        intersect=set(range(cons[first_cl]["Start"],cons[first_cl]["Stop"])).intersection(set(range(cons[second_cl]["Start"],cons[second_cl]["Stop"])))
        if only_with_common_snip==False and len(commonSNP)==0 and len(intersect)>0:
            d=0
        else:
            #for snp in SNP_pos:
            for snp in commonSNP:
                try:
                    b1=cons[first_cl][snp]
                    b2=cons[second_cl][snp]
                    if b1 != b2 and len(b1)!=0 and  len(b2)!=0:
                        if d==-1:
                            d=0
                        d=d+1
                    elif b1 == b2:
                        if d==-1:
                            d=0
                        d=d
                except:
                    continue
                if d>=1:
                    d=1
                    break

                else:
                    continue
    except(IndexError):pass
    return (d)



def build_adj_matrix_clusters (cons, SNP_pos,cl,only_with_common_snip=True):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
    #clusters=sorted(set(cl['Cluster'].values))
    Y=[]
    X=[]
    Z=[]
    sort=[]
    for k,v in cons.items():
        X.append(k)
        Y.append(int(v["Start"]))
        Z.append(int(v["Stop"]))
        sort.append([k,int(v["Start"]),int(v["Stop"])])


    from more_itertools import sort_together
    #sorted_by_pos=sort_together([Y, X])[1]
    sorted_by_pos=[]
    for i in sorted(sort, key=lambda sort: [sort[1], sort[2]]):
        sorted_by_pos.append(i[0])
    clusters=sorted(set(sorted_by_pos) & set(clusters), key=sorted_by_pos.index)
    #print(clusters)
    m = pd.DataFrame(-1, index=clusters, columns=clusters)
    for i in range(0,m.shape[1]):
        first_cl=m.index[i]
        for k in range(i+1,m.shape[1]):
            second_cl=m.index[k]
            #try:
            if m[second_cl][first_cl]==-1:
                if only_with_common_snip==True:
                    m[second_cl][first_cl] = distance_clusters(first_cl,second_cl, cons,SNP_pos)
                else:
                    m[second_cl][first_cl] = distance_clusters(first_cl, second_cl, cons, SNP_pos,False)
            #except: (KeyError)
    return (m)


def join_clusters(cons, SNP_pos, cl,R, edge, only_with_common_snip=True):
    if only_with_common_snip==False:
        M = build_adj_matrix_clusters(cons, SNP_pos, cl, False)
    else:
        M = build_adj_matrix_clusters(cons, SNP_pos, cl)


    M=change_w(M,R)
    G_vis = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    to_remove = []
    G_vis_temp = nx.nx_agraph.to_agraph(G_vis)
    G_vis_temp.layout(prog="neato")
    G_vis_temp.draw("output/graphs/cluster_GV_graph_%s_beforeremove.png" % (edge))

    path_remove=[]
    for node in G_vis.nodes():
        neighbors = nx.all_neighbors(G_vis, node)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G_vis, node, neighbor):
                if len(n_path) == 3:
                    path_remove.append(n_path)

    for n_path in path_remove:
        try:
            G_vis.remove_edge(n_path[0], n_path[2])
        except:
            continue

    lis = list(nx.topological_sort(nx.line_graph(G_vis)))
    first = []
    last = []

    for i in lis:
        first.append(i[0])
        last.append(i[1])

    for i in lis:
        if first.count(i[0]) > 1 or last.count(i[1]) > 1:
            to_remove.append(i)


    G_vis.remove_edges_from(ebunch=to_remove)


    G_vis = nx.nx_agraph.to_agraph(G_vis)
    G_vis.layout(prog="neato")
    G_vis.draw("output/graphs/cluster_GV_graph_%s.png" % (edge))


    G = nx.from_pandas_adjacency(M)
    for n_path in path_remove:
        try:
            G.remove_edge(n_path[0], n_path[2])
        except :
            continue
    G.remove_edges_from(ebunch=to_remove)
    groups = list(nx.connected_components(G))


    #print(groups)
    for group in groups:
        if len(group) > 0:
            for i in range(1, len(group)):
                cl.loc[cl['Cluster'] == list(group)[i], 'Cluster'] = list(group)[0]
    return (cl)

def postprocess (bam,cl,SNP_pos, data, edge, R, I):
    cons=build_data_cons(cl, SNP_pos, data)
    for key, val in cons.copy().items():
        if val["Strange"]==1:
            cluster=key
            clSNP=val["clSNP"]
            child_clusters=[]
            split_cluster(cl, cluster, data, clSNP, bam, edge,child_clusters, R, I)
            for child in child_clusters:
                cluster=child
                cluster_consensuns(cl,cluster,SNP_pos, data, cons)

    if len(cl.loc[cl['Cluster'] == 1000000]['ReadName'].values) != 0:
        cluster_consensuns(cl, 1000000, SNP_pos, data, cons)
        cluster = 1000000
        val=cons[cluster]
        clSNP = SNP_pos
        child_clusters = []
        #child_clusters = []
        split_cluster(cl, cluster, data, clSNP, bam, edge, child_clusters, R, I)
        for child in child_clusters:
            cluster = child
            cluster_consensuns(cl, cluster, SNP_pos, data, cons)


    cl.to_csv("output/clusters/clusters_before_joining_%s_%s_%s.csv" % (edge, I, 0.1))
    cl=join_clusters(cons, SNP_pos, cl,R, edge)
    #cl = join_clusters(cons, SNP_pos, cl, R, edge, False)
    return(cl)

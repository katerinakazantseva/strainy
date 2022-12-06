import csv
import pysam
import pandas as pd
from collections import Counter
import sys,os,subprocess
from params import *
import networkx as nx
import matplotlib.pyplot as plt
from build_adj_matrix import *
from build_data  import *
import pygraphviz as gv
import pylab
from community_detection import find_communities
from more_itertools import sort_together



def split_cluster(cl,cluster, data,clSNP, bam, edge, child_clusters, R, I,only_with_common_snip=True):
    #print("Split cluster  "+str(cluster))
    reads=sorted(set(cl.loc[cl['Cluster'] == cluster]['ReadName'].values))
    if cluster==unclustered_group_N or cluster==unclustered_group_N2  or only_with_common_snip==False: #NA cluster
    #if cluster == unclustered_group_N and len(clSNP)==0:  # NA cluster
        #print("Build na matrix")
        m = build_adj_matrix(cl[cl['Cluster'] == cluster], data, clSNP, I, bam, edge, R, only_with_common_snip=False)
    else:
        m=build_adj_matrix(cl[cl['Cluster'] == cluster], data, clSNP, I, bam,edge,R)
    m = remove_edges(m, 1)
    m.columns=range(0,len(cl[cl['Cluster'] == cluster]['ReadName']))
    m.index=range(0,len(cl[cl['Cluster'] == cluster]['ReadName']))
    m = change_w(m, R)
    G_sub = nx.from_pandas_adjacency(m)



    cluster_membership = find_communities(G_sub)
    clN=0
    uncl=0
    reads = cl[cl['Cluster'] == cluster]['ReadName'].values
    #print(cluster_membership)
    for value in set(cluster_membership.values()):
        group = [k for k, v in cluster_membership.items() if v == value]
        if len(group) > min_cluster_size:
            clN = clN + 1
            #print("NEW cluster "+str(cluster+split_id+clN))
            child_clusters.append(cluster+split_id+clN)
            for i in group:
                cl.loc[cl['ReadName'] == reads[i], "Cluster"] = cluster+split_id+clN

        else:
            uncl = uncl + 1
            for i in group:
                cl.loc[cl['ReadName'] == reads[i], "Cluster"] = unclustered_group_N2

                '''if cluster==unclustered_group_N:
                    cl.loc[cl['ReadName'] == reads[i], "Cluster"] = 'NA'
                    #cl.loc[cl['ReadName'] == reads[i], "Cluster"] = unclustered_group_N
                    #child_clusters.append(cluster + split_id + clN)
                else:
                    cl.loc[cl['ReadName'] == reads[i], "Cluster"] = unclustered_group_N
                    #child_clusters.append(cluster + split_id + clN)
                    #cl.loc[cl['ReadName'] == reads[i], "Cluster"] = cluster + split_id
                '''


    #print("New clusters  "+str(clN))
    #print("Unclustered"+str(uncl))







def build_adj_matrix_clusters (cons, SNP_pos,cl,only_with_common_snip=True):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
    try:
        clusters.remove(0)
    except:
        pass
    Y=[]
    X=[]
    Z=[]
    sort=[]

    for k,v in cons.items():
        X.append(k)
        Y.append(int(v["Start"]))
        Z.append(int(v["Stop"]))
        sort.append([k,int(v["Start"]),int(v["Stop"])])
    sorted_by_pos=[]

    for i in sorted(sort, key=lambda sort: [sort[1], sort[2]]):
        sorted_by_pos.append(i[0])
    clusters=sorted(set(sorted_by_pos) & set(clusters), key=sorted_by_pos.index)
    m = pd.DataFrame(-1, index=clusters, columns=clusters)

    for i in range(0,m.shape[1]):
        first_cl=m.index[i]
        for k in range(i+1,m.shape[1]):
            second_cl=m.index[k]

            if m[second_cl][first_cl]==-1:
                if only_with_common_snip==True:
                    m[second_cl][first_cl] = distance_clusters(first_cl,second_cl, cons,SNP_pos)
                else:
                    m[second_cl][first_cl] = distance_clusters(first_cl, second_cl, cons, SNP_pos,False)
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
    G_vis_before = nx.nx_agraph.to_agraph(G_vis)
    G_vis_before.layout(prog="neato")
    G_vis_before.draw("%s/graphs/cluster_GV_graph_before_remove_%s.png" % (output,edge))

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
    G_vis.draw("%s/graphs/cluster_GV_graph_%s.png" % (output,edge))
    G = nx.from_pandas_adjacency(M)
    for n_path in path_remove:
        try:
            G.remove_edge(n_path[0], n_path[2])
        except :
            continue
    G.remove_edges_from(ebunch=to_remove)
    groups = list(nx.connected_components(G))
    for group in groups:
        if len(group) > 0:
            for i in range(1, len(group)):
                #print(type(list(group)[0]))
                cl.loc[cl['Cluster'] == list(group)[i], 'Cluster'] = int(list(group)[0])
    return (cl)

def postprocess (bam,cl,SNP_pos, data, edge, R, I):

    cons=build_data_cons(cl, SNP_pos, data,edge)
    for key, val in cons.copy().items():
        if val["Strange"]==1 and key!=unclustered_group_N:
            cluster=key
            clSNP=val["clSNP"]
            child_clusters=[]
            split_cluster(cl, cluster, data, clSNP, bam, edge,child_clusters, R, I)
            for child in set(child_clusters):
                cluster_consensuns(cl,child,SNP_pos, data, cons,edge)
                if cons[child]["Strange"]==1:
                    split_cluster(cl, child, data, clSNP, bam, edge, child_clusters, R, I)



    cluster=unclustered_group_N
    cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge)
    child_clusters = []
    clSNP = cons[cluster]["clSNP"]
    split_cluster(cl, cluster, data, clSNP, bam, edge, child_clusters, R, I)

    for child in set(child_clusters):
        cluster_consensuns(cl, child, SNP_pos, data, cons, edge)
        if cons[child]["Strange"] == 1:
            split_cluster(cl, child, data, clSNP, bam, edge, child_clusters, R, I)

    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))

    print("Split2")
    for cluster in clusters:
        try:
            if cons[cluster]["Strange2"]==1 and cluster!=unclustered_group_N:
                clSNP=cons[cluster]["clSNP2"]
                #cluster = key
                child_clusters=[]
                split_cluster(cl, cluster, data, clSNP, bam, edge,child_clusters, R, I)
                for child in set(child_clusters):
                    cluster_consensuns(cl,child,SNP_pos, data, cons,edge)
                    if cons[child]["Strange2"]==1:
                        split_cluster(cl, child, data, clSNP, bam, edge, child_clusters, R, I)
        except(KeyError):
            continue



    if len(cl.loc[cl['Cluster'] == unclustered_group_N2]['ReadName'].values) != 0:
        cluster_consensuns(cl, unclustered_group_N2, SNP_pos, data, cons,edge)
        cluster = unclustered_group_N2
        val=cons[cluster]
        clSNP = SNP_pos
        #clSNP = val["clSNP2"]
        child_clusters = []
        split_cluster(cl, cluster, data, clSNP, bam, edge, child_clusters, R, I,only_with_common_snip = False)
        for child in set(child_clusters):
            cluster = child
            cluster_consensuns(cl, cluster, SNP_pos, data, cons,edge)
            split_cluster(cl, cluster, data, clSNP, bam, edge, child_clusters, R, I,only_with_common_snip = False)

    #clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
    #clusters = sorted(set(cl.loc[cl['Cluster'] != unclustered_group_N2]['Cluster'].values))
    cl = cl[cl['Cluster'] != unclustered_group_N2]

    counts = cl['Cluster'].value_counts(dropna=False)
    cl = cl[~cl['Cluster'].isin(counts[counts < 6].index)]  #change for cov*01.

    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))


    cl.to_csv("%s/clusters/clusters_before_joining_%s_%s_%s.csv" % (output,edge, I, 0.1))

    cons = build_data_cons(cl, SNP_pos, data, edge)
    cl=join_clusters(cons, SNP_pos, cl,R, edge)
    cons = build_data_cons(cl, SNP_pos, data, edge)
    cl=join_clusters(cons, SNP_pos, cl, R, edge, False)
    counts = cl['Cluster'].value_counts(dropna=False)
    cl = cl[~cl['Cluster'].isin(counts[counts < 6].index)] #change for cov*01.
    return(cl)

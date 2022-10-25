import networkx as nx
from community_detection import find_communities

from build_adj_matrix import *
from build_data import *


def split_cluster(cl, cluster, data, clSNP, bam, edge, child_clusters, R, I):
    print("Strange cluster detected")
    if cluster==1000000:
        print("Build na matrix")
        m = build_adj_matrix(cl[cl['Cluster'] == cluster], data, clSNP, I, bam, edge, R, only_with_common_snip=False)
    else:
        m = build_adj_matrix(cl[cl['Cluster'] == cluster], data, clSNP, I, bam, edge, R)
    m = remove_edges(m, 1)
    m.columns=range(0,len(cl[cl['Cluster'] == cluster]['ReadName']))
    m.index=range(0,len(cl[cl['Cluster'] == cluster]['ReadName']))
    m = change_w(m, R)
    G = nx.from_pandas_adjacency(m)
    cluster_membership = find_communities(G)
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


def build_adj_matrix_clusters (cons, cl, edge, flye_consensus, only_with_common_snip=True):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
    sort=[]
    for k,v in cons.items():
        sort.append([k, int(v["Start"]), int(v["Stop"])])
    sorted_by_pos=[]
    for i in sorted(sort, key=lambda sort: [sort[1], sort[2]]):
        sorted_by_pos.append(i[0])
    clusters = sorted(set(sorted_by_pos) & set(clusters), key=sorted_by_pos.index)
    m = pd.DataFrame(-1, index=clusters, columns=clusters)
    for i in range(0, m.shape[1]):
        first_cl = m.index[i]
        for k in range(i+1,m.shape[1]):
            second_cl = m.index[k]
            if m[second_cl][first_cl] == -1:
                m[second_cl][first_cl] = flye_consensus.cluster_distance_via_alignment(first_cl, second_cl, cl, edge)
    return m


def join_clusters(cons, cl, R, edge, consensus, only_with_common_snip=True):
    if only_with_common_snip==False:
        M = build_adj_matrix_clusters(cons, cl, edge, consensus, False)
    else:
        M = build_adj_matrix_clusters(cons, cl, edge, consensus)

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
    for group in groups:
        if len(group) > 0:
            for i in range(1, len(group)):
                cl.loc[cl['Cluster'] == list(group)[i], 'Cluster'] = list(group)[0]
    return cl


def postprocess (bam, cl, SNP_pos, data, edge, R, I, flye_consensus):
    cons=build_data_cons(cl, SNP_pos, data)
    for key, val in cons.copy().items():
        if val["Strange"]==1:
            cluster=key
            clSNP=val["clSNP"]
            child_clusters=[]
            split_cluster(cl, cluster, data, clSNP, bam, edge, child_clusters, R, I)
            for child in child_clusters:
                cluster=child
                cluster_consensuns(cl,cluster,SNP_pos, data, cons)

    if len(cl.loc[cl['Cluster'] == 1000000]['ReadName'].values) != 0:
        cluster_consensuns(cl, 1000000, SNP_pos, data, cons)
        cluster = 1000000
        clSNP = SNP_pos
        child_clusters = []
        split_cluster(cl, cluster, data, clSNP, bam, edge, child_clusters, R, I)
        for child in child_clusters:
            cluster = child
            cluster_consensuns(cl, cluster, SNP_pos, data, cons)

    cl.to_csv("output/clusters/clusters_before_joining_%s_%s_%s.csv" % (edge, I, 0.1))
    cl = join_clusters(cons, cl, R, edge, flye_consensus)
    return cl

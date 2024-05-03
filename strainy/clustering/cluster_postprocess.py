import networkx as nx
import logging
import pandas as pd

from strainy.clustering.community_detection import find_communities
import strainy.clustering.build_adj_matrix as matrix
import strainy.clustering.build_data as build_data
import strainy.gfa_operations.gfa_ops as gfa_ops
from strainy.params import *

logger = logging.getLogger()


def split_cluster(cl,cluster, data,cons,clSNP, bam, edge, R, I,only_with_common_snip=True):
    #logging.debug("Split cluster: " + str(cluster)+ " "+ str(only_with_common_snip))
    child_clusters = []
    reads = sorted(set(cl.loc[cl["Cluster"] == cluster,"ReadName"].values))
    if cluster == UNCLUSTERED_GROUP_N or cluster==UNCLUSTERED_GROUP_N2  or only_with_common_snip==False: #NA cluster
        m = matrix.build_adj_matrix(cl[cl["Cluster"] == cluster], data, clSNP, I, bam, edge, R, only_with_common_snip=False)
    else:
        m = matrix.build_adj_matrix(cl[cl["Cluster"] == cluster], data, clSNP, I, bam,edge,R)
    m = matrix.remove_edges(m, 1)
    m.columns = range(0,len(cl[cl["Cluster"] == cluster]["ReadName"]))
    m.index = range(0,len(cl[cl["Cluster"] == cluster]["ReadName"]))
    m = matrix.change_w(m, R)
    G_sub = gfa_ops.from_pandas_adjacency_notinplace(m)
    cl_exist = sorted(set(cl.loc[cl["Cluster"] != "NA","Cluster"].values))+list(cons.keys())
    cluster_membership = find_communities(G_sub)
    clN = 0
    uncl = 0
    reads = cl[cl["Cluster"] == cluster]["ReadName"].values
    new_cl_id_na = cluster + SPLIT_ID
    while new_cl_id_na in cl_exist:
        new_cl_id_na = new_cl_id_na + 1
    if len(set(cluster_membership.values())) > 0:
        for value in set(cluster_membership.values()):
            group = [k for k, v in cluster_membership.items() if v == value]
            if len(group) > min_cluster_size:
                clN = clN + 1
                new_cl_id = new_cl_id_na+clN
                while new_cl_id in cl_exist:
                    new_cl_id = new_cl_id+1
                    child_clusters.append(new_cl_id)
                cl_exist.append(new_cl_id)
                for i in group:
                     mask = cl["ReadName"] == str(reads[i])
                     cl.loc[mask, "Cluster"] = new_cl_id
            else:
                uncl = uncl + 1
                for i in group:
                    mask = cl["ReadName"] == str(reads[i])
                if only_with_common_snip == True and cluster!=UNCLUSTERED_GROUP_N: #change it for parameter process NA or not
                    cl.loc[mask, "Cluster"] = new_cl_id_na
                    child_clusters.append(new_cl_id_na)
                else:
                    cl.loc[mask, "Cluster"] = UNCLUSTERED_GROUP_N
    return [new_cl_id_na, clN]


def build_adj_matrix_clusters(edge,cons,cl,flye_consensus, only_with_common_snip=True, set_slusters=None):
    if set_slusters == None:
        clusters = sorted(set(cl.loc[cl["Cluster"] != "NA","Cluster"].values))
    else:
        clusters = set_slusters
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
        Z.append(int(v["End"]))
        sort.append([k,int(v["Start"]),int(v["End"])])
    sorted_by_pos = []

    for i in sorted(sort, key = lambda sort: [sort[2], sort[1]]):
        sorted_by_pos.append(i[0])
    clusters = sorted(set(sorted_by_pos) & set(clusters), key = sorted_by_pos.index)
    m = pd.DataFrame(-1.0, index = clusters, columns = clusters)

    for i in range(0, m.shape[1]):
        first_cl = m.index[i]
        for k in range(i + 1, m.shape[1]):
            second_cl = m.index[k]

            if m[second_cl][first_cl] == -1:
                m.loc[first_cl, second_cl] = matrix.distance_clusters(edge, first_cl, second_cl, cons, cl,flye_consensus, only_with_common_snip)
                #m[second_cl][first_cl] = matrix.distance_clusters(edge, first_cl, second_cl, cons, cl,flye_consensus, only_with_common_snip)
    return m



def join_clusters(cons, cl, Rcl, edge, consensus, only_with_common_snip=True,set_clusters=None, only_nested=False, transitive=False):
    MAX_VIS_SIZE = 500
    CUT_OFF=3

    if only_with_common_snip == False:
        if set_clusters == None:
            M = build_adj_matrix_clusters(edge,cons, cl,consensus, False)
        else:
            M = build_adj_matrix_clusters(edge, cons, cl, consensus, False,set_clusters)
    else:
        if set_clusters == None:
            M = build_adj_matrix_clusters(edge,cons, cl,consensus, True)
        else:
            M = build_adj_matrix_clusters(edge, cons, cl, consensus, True, set_clusters)
    M = matrix.change_w(M,Rcl)

    try:
        G_vis = gfa_ops.from_pandas_adjacency_notinplace(M, create_using = nx.DiGraph)
    except nx.NetworkXUnfeasible:
        M.index = M.index.astype(str)
        M.columns = M.columns.astype(str)
        G_vis = gfa_ops.from_pandas_adjacency_notinplace(M, create_using = nx.DiGraph)

    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))

    if StRainyArgs().debug and max(G_vis.number_of_nodes(), G_vis.number_of_edges()) < MAX_VIS_SIZE:
        G_vis_before = nx.nx_agraph.to_agraph(G_vis)
        G_vis_before.layout(prog = "dot")
        G_vis_before.draw("%s/graphs/linear_phase_%s.png" % (StRainyArgs().output_intermediate, edge))
        CUT_OFF=5

    path_remove = []
    for node in G_vis.nodes():
        neighbors = nx.all_neighbors(G_vis, node)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G_vis, node, neighbor, cutoff=CUT_OFF):
                if len(n_path) >2:
                    path_remove.append(n_path)

    for n_path in path_remove:
        try:
            G_vis.remove_edge(n_path[0], n_path[len(n_path)-1])
        except:
            continue

    lis = list(nx.topological_sort(nx.line_graph(G_vis)))
    first = []
    last = []

    for i in lis:
        first.append(i[0])
        last.append(i[1])

    to_remove = []
    for i in lis:
        if first.count(i[0]) > 1 or last.count(i[1]) > 1:
            to_remove.append(i)
    G_vis.remove_edges_from(ebunch = to_remove)

    if StRainyArgs().debug and max(G_vis.number_of_nodes(), G_vis.number_of_edges()) < MAX_VIS_SIZE:
        G_vis = nx.nx_agraph.to_agraph(G_vis)
        G_vis.layout(prog="dot")
        G_vis.draw("%s/graphs/linear_phase_simplified_%s.png" % (StRainyArgs().output_intermediate, edge))

    G = gfa_ops.from_pandas_adjacency_notinplace(M)


    if transitive==True:
        graph_str = str(nx.nx_agraph.to_agraph(G))
        not_visited=list(G.nodes).copy()
        while not_visited:
             
            node=not_visited[0]
            neighbors = nx.all_neighbors(G, node)
            n_list=list(neighbors)
            for nei in n_list:
                n_list2=list(nx.all_neighbors(G, nei))
                n_list.remove(nei)
                n_list2.remove(node)
                if (n_list==n_list2) is True:
                    сlusters = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))
                    new_cl_id = int(nei)+UNCLUSTERED_GROUP_N
                    while new_cl_id in сlusters:
                        new_cl_id = new_cl_id+1
                    cl.loc[cl["Cluster"] == int(node), "Cluster"] =new_cl_id
                    cl.loc[cl["Cluster"] == int(nei), "Cluster"] =new_cl_id

                    G.remove_node(node)
                    G=nx.relabel_nodes(G,{nei:new_cl_id})

                    not_visited.append(new_cl_id)
                    try:
                        not_visited.remove(nei)
                        not_visited.remove(node)
                    except:
                        pass
                    break
            try:
                not_visited.remove(node)
            except:
                 pass


    else:
        for n_path in path_remove:
            try:
                G.remove_edge(n_path[0], n_path[len(n_path)-1])
            except :
                continue
        G.remove_edges_from(ebunch = to_remove)

        nested = {}
        nodes = list(G.nodes())




        for node in nodes:
            try:
                neighbors = nx.all_neighbors(G, node)
                for neighbor in list(neighbors):
                    if cons[int(node)]["Start"] < cons[int(neighbor)]["Start"] and cons[int(node)]["End"] > cons[int(neighbor)]["End"]:
                        try:
                            G.remove_edge(node, neighbor)
                            logger.debug("REMOVE NESTED" + str(neighbor)+" :"+str(node))
                            if len(nx.all_neighbors(G, neighbor)) == 1:
                                try:
                                    nested[neighbor] = nested[neighbor].append(node)
                                except:
                                    nodes = [node]
                                    nested[neighbor] = nodes
                        except:
                            continue
            except:
                continue

        groups = list(nx.connected_components(G))
        if only_nested == True:
            for k,v in nested.items():
                if len(v) == 1:
                    cl.loc[cl["Cluster"] == int(k), "Cluster"] = int(v[0])

        else:
            for group in groups:
                if len(group) > 1:
                    сlusters = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))
                    new_cl_id = int(list(group)[0])+UNCLUSTERED_GROUP_N
                    while new_cl_id in сlusters:
                        new_cl_id = new_cl_id+1
                    for i in range(0, len(group)):
                        cl.loc[cl["Cluster"] == int(list(group)[i]), "Cluster"] = new_cl_id
    return cl


def split_all(cl, cluster, data, cons,bam, edge, R, I, SNP_pos,reference_seq,type):
    if type=="unclustered":
        factor="Strange"
        snp_set="clSNP"
    if type=="lowheterozygosity":
        factor="Strange2"
        snp_set="clSNP2"


    if cons[cluster][factor] == 1:
        clSNP = cons[cluster][snp_set]
        res = split_cluster(cl, cluster, data,cons, clSNP, bam, edge, R, I)
        new_cl_id_na=res[0]
        clN =res[1]
        build_data.cluster_consensuns(cl, new_cl_id_na, SNP_pos, data, cons, edge, reference_seq)

        if clN != 0: #if clN==0 we dont need split NA cluster
            split_cluster(cl, new_cl_id_na, data, cons,cons[new_cl_id_na][snp_set], bam, edge, R, I, False)
        clusters = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))

        if clN == 1: #STOP LOOP IF EXIST
            build_data.cluster_consensuns(cl, new_cl_id_na + clN, SNP_pos, data, cons, edge, reference_seq)

        for cluster in clusters:
            if cluster not in cons:
                build_data.cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge, reference_seq)
                split_all(cl, cluster, data, cons,bam, edge, R, I, SNP_pos,reference_seq,"unclustered")






def postprocess(bam, cl, SNP_pos, data, edge, R,Rcl, I, flye_consensus,mean_edge_cov):
    reference_seq = build_data.read_fasta_seq(StRainyArgs().fa, edge)
    cons = build_data.build_data_cons(cl, SNP_pos, data, edge, reference_seq)
    if StRainyArgs().debug:
        cl.to_csv("%s/clusters/%s_1.csv" % (StRainyArgs().output_intermediate, edge))
    clusters = sorted(set(cl.loc[cl["Cluster"] != "NA","Cluster"].values))



    cl.loc[cl["Cluster"] == "NA", "Cluster"] = UNCLUSTERED_GROUP_N
    build_data.cluster_consensuns(cl, UNCLUSTERED_GROUP_N, SNP_pos, data, cons, edge, reference_seq)
    clSNP = cons[UNCLUSTERED_GROUP_N]["clSNP2"]
    splitna = split_cluster(cl, UNCLUSTERED_GROUP_N, data, cons, clSNP, bam, edge, R, I,False)

    #Remove unclustered reads after splitting NA cluster
    splitna[0]
    cl = cl[cl["Cluster"] != splitna[0]]

    cl = cl[cl["Cluster"] != UNCLUSTERED_GROUP_N]
    clusters = sorted(set(cl.loc[cl["Cluster"] != splitna[0], "Cluster"].values))
    clusters = sorted(set(cl.loc[cl["Cluster"] != UNCLUSTERED_GROUP_N, "Cluster"].values))

    build_data.cluster_consensuns(cl, UNCLUSTERED_GROUP_N, SNP_pos, data, cons, edge, reference_seq)
    counts = cl["Cluster"].value_counts(dropna = False)

    for cluster in clusters:
        if cluster not in cons:
            build_data.cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge, reference_seq)

    cl = join_clusters(cons, cl, Rcl, edge, flye_consensus)
    cons = build_data.build_data_cons(cl, SNP_pos, data, edge, reference_seq)

    clusters = sorted(set(cl.loc[cl["Cluster"] != "NA","Cluster"].values))
    prev_clusters = clusters
    for cluster in clusters:
        split_all(cl, cluster, data, cons,bam, edge, R, I, SNP_pos,reference_seq,"unclustered")
        clusters = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))
        new_clusters = list(set(clusters) - set(prev_clusters))
        prev_clusters = clusters
        cl= join_clusters(cons, cl, Rcl, edge, flye_consensus, False,new_clusters,only_nested=False)
        clusters = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))

        for cluster in clusters:
            if cluster not in cons:
                build_data.cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge, reference_seq)
    clusters = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))

    logging.info("Split stage2: Break regions of low heterozygosity")
    for cluster in clusters:
        split_all(cl, cluster, data, cons,bam, edge, R, I, SNP_pos,reference_seq,"lowheterozygosity")




    cl.loc[cl["Cluster"] == "NA", "Cluster"] = UNCLUSTERED_GROUP_N
    build_data.cluster_consensuns(cl, UNCLUSTERED_GROUP_N, SNP_pos, data, cons, edge, reference_seq)
    clSNP = cons[UNCLUSTERED_GROUP_N]["clSNP2"]
    splitna = split_cluster(cl, UNCLUSTERED_GROUP_N, data, cons, clSNP, bam, edge, R, I,False)
    
    #Remove unclustered reads after splitting NA cluster
    splitna[0]
    cl = cl[cl["Cluster"] != splitna[0]]

    cl = cl[cl["Cluster"] != UNCLUSTERED_GROUP_N]
    clusters = sorted(set(cl.loc[cl["Cluster"] != splitna[0], "Cluster"].values))
    clusters = sorted(set(cl.loc[cl["Cluster"] != UNCLUSTERED_GROUP_N, "Cluster"].values))

    cl=update_cluster_set(cl, cluster, SNP_pos, data, cons, edge, reference_seq,mean_edge_cov)

    cl = join_clusters(cons, cl, Rcl, edge, flye_consensus)
    cl=update_cluster_set(cl, cluster, SNP_pos, data, cons, edge, reference_seq,mean_edge_cov)
    cl = join_clusters(cons, cl, Rcl, edge, flye_consensus, transitive=True)
    cl=update_cluster_set(cl, cluster, SNP_pos, data, cons, edge, reference_seq,mean_edge_cov)
    cl = join_clusters(cons, cl, Rcl, edge, flye_consensus)
    cl=update_cluster_set(cl, cluster, SNP_pos, data, cons, edge, reference_seq,mean_edge_cov,fraction=0.05)
    cl = join_clusters(cons, cl, Rcl, edge, flye_consensus, only_with_common_snip=False,only_nested=True)
    counts = cl["Cluster"].value_counts(dropna = False)
    cl = cl[~cl["Cluster"].isin(counts[counts < 6].index)]  #TODO change for cov*01.
    cl=update_cluster_set(cl, cluster, SNP_pos, data, cons, edge, reference_seq,mean_edge_cov,fraction=0.05)
    return cl



def update_cluster_set(cl, cluster, SNP_pos, data, cons, edge, reference_seq,mean_edge_cov,fraction=0.01):
    #Update consensus and remove small clusters (less 5% of unitig coverage)
    clusters = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))
    for cluster in clusters:
        if cluster not in cons:
            build_data.cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge, reference_seq)

    for cluster in clusters:
        if cons[cluster]['Cov']<mean_edge_cov*fraction:
            cl = cl[cl["Cluster"] !=cluster]


    return cl

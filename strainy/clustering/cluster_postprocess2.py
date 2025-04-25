import networkx as nx
import logging
import pandas as pd
from strainy.clustering.community_detection import find_communities
from strainy.clustering import build_adj_matrix as matrix
from strainy.clustering import build_data
from strainy.graph_operations import gfa_ops
from strainy.params import *
from functools import partial
from multiprocessing import Pool
from tqdm import tqdm



logger = logging.getLogger()

'''
def overlap(cons, cl, Rcl, data, snp_pos, edges, consensus):
    M = build_adj_matrix_clusters_whole(edges, cons, cl, consensus, data,snp_pos, False)
    try:
        G_vis = gfa_ops.from_pandas_adjacency_notinplace(M, create_using = nx.DiGraph)
    except nx.NetworkXUnfeasible:
        M.index = M.index.astype(str)
        M.columns = M.columns.astype(str)
        G_vis = gfa_ops.from_pandas_adjacency_notinplace(M, create_using = nx.DiGraph)
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    MAX_VIS_SIZE = 500
    CUT_OFF=3
    if max(G_vis.number_of_nodes(), G_vis.number_of_edges()) < MAX_VIS_SIZE:
        G_vis_before = nx.nx_agraph.to_agraph(G_vis)
        G_vis_before.layout(prog = "dot")
        G_vis_before.draw(f"{StRainyArgs().output_intermediate}/graphs/linear_phase_{edge}.png")
        CUT_OFF=5

'''

def pairs(edges,cl,flye_consensus,data, snp_pos,reversed_edges,
                              only_with_common_snip=True):

    clusters = sorted(set(cl.loc[cl["Cluster"] != "NA","Cluster"].values))
    try:
        clusters.remove(0)
    except ValueError:
        pass

    m = pd.DataFrame(-1.0, index = clusters, columns = clusters)
    pairs=[]
    for i in range(0, m.shape[1]):
        first_cl = m.index[i]
        for k in range(i + 1, m.shape[1]):
            second_cl = m.index[k]
            if m[second_cl][first_cl] == -1:
                cluster_d=matrix.distance_clusters2\
                    (edges, first_cl, second_cl, data, cl,flye_consensus,snp_pos, only_with_common_snip)
                m.loc[first_cl, second_cl] = cluster_d
                if cluster_d==0:
                    pairs.append(matrix.orient(edges, first_cl,second_cl, cl, reversed_edges))
    return pairs



def pairs_parallel(edges, cl, flye_consensus, data, snp_pos, reversed_edges,
                   output_intermediate, only_with_common_snip=True, n_processes=4):
    clusters = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))
    try:
        clusters.remove(0)
    except ValueError:
        pass

    m = pd.DataFrame(-1.0, index=clusters, columns=clusters)
    cluster_pairs = []

    for i in range(0, m.shape[1]):
        first_cl = m.index[i]
        for k in range(i + 1, m.shape[1]):
            second_cl = m.index[k]
            if m[second_cl][first_cl] == -1:
                cluster_d = matrix.distance_clusters2(
                    edges, first_cl, second_cl, data, cl, flye_consensus, snp_pos, only_with_common_snip
                )
                m.loc[first_cl, second_cl] = cluster_d
                if cluster_d == 0:
                    cluster_pairs.append((first_cl, second_cl))

    chunk_size = len(cluster_pairs) // n_processes + 1
    chunks = [cluster_pairs[i:i + chunk_size] for i in range(0, len(cluster_pairs), chunk_size)]

    with Pool(processes=n_processes) as pool:
        func = partial(process_pair, edges=edges, cl=cl, flye_consensus=flye_consensus, data=data,
                       snp_pos=snp_pos, reversed_edges=reversed_edges, output_intermediate=output_intermediate)
        results = list(tqdm(pool.imap(func, chunks), total=len(chunks), desc="Cluster distances"))

    return results


def process_pair(chunk, edges, cl, flye_consensus, data, snp_pos, reversed_edges, output_intermediate):
    results = []
    for first_cl, second_cl in chunk:
        flye_dir = f"{output_intermediate}/flye_outputs/asm/{first_cl}/asm"
        cluster_d = matrix.distance_clusters2(
            edges, first_cl, second_cl, data, cl, flye_consensus, snp_pos, only_with_common_snip=True
        )
        results.append((first_cl, second_cl, cluster_d))
    return results
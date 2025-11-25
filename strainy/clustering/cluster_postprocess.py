import networkx as nx
import logging
import pandas as pd
from strainy.clustering.community_detection import find_communities
from strainy.clustering import build_adj_matrix as matrix
from strainy.clustering import build_data
from strainy.graph_operations import gfa_ops
from strainy.params import *
from tqdm import tqdm
logger = logging.getLogger()


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
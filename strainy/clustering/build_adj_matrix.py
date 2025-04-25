import logging
import pandas as pd
from scipy.spatial.distance import cdist
from strainy.params import *
from strainy.clustering import build_data
import pandas as pd
import numpy as np
import tempfile
import os
logger = logging.getLogger()
import networkx as nx
import networkit as nk
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from functools import partial
#pd.options.mode.chained_assignment = None
from collections import Counter


class DistanceWrapper2():
    # Wrapper for calling cdist with custom distance function
    def __init__(self, edges, data, snp_pos, R, only_with_common_snip):
        self.data = data
        self.snp_pos = snp_pos
        self.R = R
        self.edges=edges
        self.only_with_common_snip = only_with_common_snip

    def distance_wrapper2(self, first_read, second_read):
        """
        return distance2(first_read[0],
                        second_read[0],
                        self.data,
                        self.snp_pos,
                        self.R,
                        self.edges,
                        self.only_with_common_snip)
        """
        return distance2(first_read,
                        second_read,
                        self.data,
                        self.snp_pos,
                        self.R,
                        self.edges,
                        self.only_with_common_snip)






def distance_clusters2(edges,first_cl,second_cl, data,cl, flye_consensus,snp_pos, only_with_common_snip=True):
    """
    Calculates the distance between two clusters based on SNP data and sequence alignment.
    This function computes a distance measure between two clusters (`first_cl` and `second_cl`) on a given
    `edge`. The distance is based on the overlap of SNPs, the number of common SNPs, and sequence alignment,
    and it provides a normalized value representing how similar or different the clusters are.
    Returns:
        float: A floating-point value representing the distance between the two clusters.
               - A value of `0` indicates high similarity.
               - A value of `1` indicates no similarity or no significant overlap.
    """
    #TODO cacl distances for all edges (first we need to check common)
    #d = -1
    intersect=0
    for edge in edges:
        bound_first = build_data.cluster_bounds(cl, first_cl, edge)
        bound_second = build_data.cluster_bounds(cl, second_cl, edge)
        try:
            intersect_edge = min(bound_first[1], bound_second[1]) - max(bound_first[0], bound_second[0])
            if intersect_edge<0:
                intersect_edge=0
            intersect=intersect+intersect_edge
        except KeyError:
            continue
    if intersect > I:
        dc=0
        for edge in edges:
            bound_first = build_data.cluster_bounds(cl, first_cl, edge)
            bound_second = build_data.cluster_bounds(cl, second_cl, edge)
            snips=[int(x) for x in snp_pos[edge] if int(x) < min(bound_first[1], bound_second[1]) and int(x) > max(bound_first[0], bound_second[0])]
            intersect_edge = min(bound_first[1], bound_second[1]) - max(bound_first[0], bound_second[0])
            if intersect_edge>0:
                dc=dc+(flye_consensus.cluster_distance_via_alignment(first_cl, second_cl, cl, edge, snips))
        d =dc/intersect
    else:
        d = 1
    return float(d)


def orient(edges,first_cl,second_cl,cl,reversed_edges):
    normal=tuple([first_cl,second_cl])
    reverse=tuple([second_cl,first_cl])
    #order=normal
    #TODO check only single edge
    for edge in edges:
        bound_first = build_data.cluster_bounds(cl, first_cl, edge)
        bound_second = build_data.cluster_bounds(cl, second_cl, edge)
        if bound_first!=[-1,-1] and bound_second!=[-1,-1]:
            if bound_first[0]<bound_second[0] and edge not in reversed_edges:
                order=normal
            elif bound_first[0]<bound_second[0] and edge in reversed_edges:
                order=reverse
            elif bound_first[0]>bound_second[0] and edge not in reversed_edges:
                order=reverse
            elif bound_first[0]>bound_second[0] and edge in reversed_edges:
                order=normal
            elif bound_first[1]<bound_second[1] and edge not in reversed_edges:
                order=normal
            elif bound_first[1]<bound_second[1] and edge in reversed_edges:
                order=reverse
            elif bound_first[1]>bound_second[1] and edge not in reversed_edges:
                order=reverse
            elif bound_first[1]>bound_second[1] and edge in reversed_edges:
                order=normal
    return order


"""
#non parallel
def build_adj_matrix2(edges, data, snp_pos, I, file, R,only_with_common_snip=True):
    m = pd.DataFrame(-1.0, index=list(data.keys()), columns=list(data.keys()))
    if only_with_common_snip==False:
        #dw = DistanceWrapper2(cl, data, snp_pos, R, only_with_common_snip)
        dw = DistanceWrapper2(edges, data, snp_pos, R, only_with_common_snip)
        result = cdist(pd.Series(list(data.keys()),name="ReadName").to_frame(), pd.Series(list(data.keys()),name="ReadName").to_frame(), dw.distance_wrapper2)

        # Set the first row and the column to -1
        try:
            result[0,:] = -1
            result[:,0] = -1
        except IndexError:
            pass

        result_df = pd.DataFrame(result,
                                 index=list(data.keys()),
                                 columns=list(data.keys()))
    else:
        dw = DistanceWrapper2(edges, data, snp_pos, R, True) #tofo change for arg only_with_common_snip
        result = cdist(pd.Series(list(data.keys()),name="ReadName").to_frame(), pd.Series(list(data.keys()),name="ReadName").to_frame(), dw.distance_wrapper2)
        result[0,:] = -1
        result_df = pd.DataFrame(result,
                        index=list(data.keys()),
                        columns=list(data.keys()))

    return result_df
"""    
    

#parallel
def compute_chunk(chunk, edges, data, snp_pos, R, only_with_common_snip):
    dw = DistanceWrapper2(edges, data, snp_pos, R, only_with_common_snip)
    result = []
    for a, b in chunk:
        d = dw.distance_wrapper2(a, b)
        result.append((a, b, d))
    return result


def build_adj_matrix2(
        edges, data, snp_pos, I, file, R,
        only_with_common_snip=True,
        use_memmap=True
):
    read_names = list(data.keys())
    m_size = len(read_names)
    name_to_idx = {name: idx for idx, name in enumerate(read_names)}
    if use_memmap:
        tmp_dir = tempfile.mkdtemp()
        matrix_file = os.path.join(tmp_dir, "adj_matrix.dat")
        matrix_np = np.memmap(matrix_file, dtype='float64', mode='w+', shape=(m_size, m_size))
        matrix_np.fill(-1.0)
    else:
        matrix_np = np.full((m_size, m_size), -1.0, dtype=np.float64)
    pairs = [(read_names[i], read_names[j]) for i in range(m_size) for j in range(i + 1, m_size)]
    #n_processes = min(8, cpu_count())
    n_processes=1
    chunk_size = len(pairs) // (n_processes * 4) + 1
    chunks = [pairs[i:i + chunk_size] for i in range(0, len(pairs), chunk_size)]

    with Pool(processes=n_processes) as pool:
        func = partial(compute_chunk, edges=edges, data=data, snp_pos=snp_pos, R=R,
                       only_with_common_snip=only_with_common_snip)
        for chunk_result in tqdm(pool.imap_unordered(func, chunks), total=len(chunks), desc="Building matrix"):
            for a, b, d in chunk_result:
                i, j = name_to_idx[a], name_to_idx[b]
                matrix_np[i, j] = d
                matrix_np[j, i] = d
    matrix_np[0, :] = -1
    matrix_np[:, 0] = -1
    matrix_df = pd.DataFrame(matrix_np, index=read_names, columns=read_names)

    if use_memmap:
        del matrix_np

    return matrix_df



def compute_edge_chunk(pairs, data, snp_pos, edges, R,
                       only_with_common_snip=True, weight_threshold=0):
    results = []
    for a, b in pairs:
        try:
            weight = distance2(a, b, data, snp_pos, R, edges, only_with_common_snip)
            #print(weight)
            if weight is not None and weight <= weight_threshold and weight!=-1.0:
                #print("add")
                results.append((a, b, 1))
        except Exception as e:
            continue
    return results


def distance2(read1, read2, data, snp_pos, R, edges, only_with_common_snip=True):
    """
      Calculates the distance between two reads based on shared SNP positions and sequence alignment.
      This function computes a distance measure between two reads (`read1` and `read2`) using their shared SNP
      positions. The distance is normalized based on the overlap of the reads and is used to evaluate their similarity.
      Returns:
          float: A distance measure between the two reads:
                 - `0` indicates high similarity.
                 - `-1.0` indicates insufficient overlap or no shared SNPs.
                 - Values between `0` and `1` indicate varying degrees of dissimilarity.
    """
    intersect=0
    commonSNP={}
    for edge in list(data[read1].keys()):
        try:
            snp_pos_edge=snp_pos[edge]
        except TypeError:
            snp_pos_edge=[]
        try:
            intersect_edge = max(min(data[read1][edge]["End"], data[read2][edge]["End"]) - max(data[read1][edge]["Start"], data[read2][edge]["Start"]),
                        0)
            intersect=intersect+intersect_edge
            firstSNPs = list(data[read1][edge].keys())
            secondSNPs = list(data[read2][edge].keys())
            keys = ('End', 'Start', 'Rclip', 'Lclip')
            firstSNPs = [key for key in firstSNPs if key not in keys]
            secondSNPs= [key for key in secondSNPs if key not in keys]
            commonSNP[edge] = sorted(set(firstSNPs).intersection(secondSNPs).intersection(snp_pos_edge))
        except KeyError:
            continue
    commonSNP_N=(sum([len(i) for i in(list(commonSNP.values()))]))


    d = -1
    #firstSNPs = list(data[read1].keys())
    #secondSNPs = list(data[read2].keys())
    #keys=('End','Start', 'Rclip', 'Lclip')
    #firstSNPs = [key for key in firstSNPs if key not in keys]
    #secondSNPs= [key for key in secondSNPs if key not in keys]
    #commonSNP = sorted(set(firstSNPs).intersection(secondSNPs).intersection(snp_pos))

    if read1 == read2:
        return 0

    if only_with_common_snip == True or (only_with_common_snip == False and commonSNP_N > 0):
        if intersect < I:
            return -1.0

        if commonSNP_N == 0:
            return -1.0
        for edge in commonSNP.keys():
            for snp in commonSNP[edge]:
                try:
                    b1 = data[read1][edge][snp]
                    b2 = data[read2][edge][snp]
                    if b1 != b2 and len(b1) != 0 and len(b2) != 0:
                        if d == -1:
                            d = 0
                        d = d + 1
                    elif b1 == b2:
                        if d == -1:
                            d = 0
                except:
                    continue

        d = d / intersect

    #TODO change it and test
    if commonSNP_N == 0 and only_with_common_snip == False:
        for edge in list(data[read1].keys()):
            try:
                intersect_edge = max(
                    min(data[read1][edge]["End"], data[read2][edge]["End"]) - max(data[read1][edge]["Start"],
                                                                                  data[read2][edge]["Start"]),
                    0)
                intersect = intersect + intersect_edge
            except KeyError:
                continue

        #intersect = max(min(data[read1]["End"], data[read2]["End"]) - max(data[read1]["Start"], data[read2]["Start"]), 0)
        if intersect > 0:
            d = 0
        else:
            d = 1
    return float(d)



def remove_edges(m, R):
    m_transformed = m
    m_transformed[m_transformed > R] = -1
    return m_transformed


def change_w(m, R):
    m_transformed = m
    m_transformed[m_transformed == 0] = -10
    m_transformed[m_transformed == -1] = 0
    m_transformed[m_transformed > R] = 0
    m_transformed[m_transformed == -10] = 0.000001
    return m_transformed


def pair_generator(read_names):
    m_size = len(read_names)
    for i in range(m_size):
        for j in range(i + 1, m_size):
            yield read_names[i], read_names[j]


def chunked_generator(generator, chunk_size):
    chunk = []
    for item in generator:
        chunk.append(item)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

def build_graph_from_data_networKit(data, snp_pos, edges, R, only_with_common_snip=True,
                                    weight_threshold=0.0, n_processes=4, chunk_size=10000):
    if n_processes is None:
        n_processes = min(cpu_count(), 16)

    read_names = list(data.keys())
    name_to_idx = {name: idx for idx, name in enumerate(read_names)}
    m_size = len(read_names)

    total_pairs = m_size * (m_size - 1) // 2
    G = nk.graph.Graph(n=m_size, weighted=True, directed=False)

    pair_gen = pair_generator(read_names)
    chunked_pairs = chunked_generator(pair_gen, chunk_size)

    func = partial(compute_edge_chunk, data=data, snp_pos=snp_pos, edges=edges, R=R,
                   only_with_common_snip=only_with_common_snip, weight_threshold=weight_threshold)

    with Pool(processes=n_processes) as pool:
        with tqdm(total=total_pairs, desc="Building graph", unit="pair") as pbar:
            for result in pool.imap_unordered(func, chunked_pairs):
                for a, b, weight in result:
                    G.addEdge(name_to_idx[a], name_to_idx[b], weight)
                pbar.update(len(result))
    return G


"""
def build_graph_from_data(
    edges, data, snp_pos, I, file, R,
    only_with_common_snip=True,
    use_memmap=True,
    memmap_filename=None,
    n_processes=4,
    weight_threshold=0.0:
    read_names = list(data.keys())
    m_size = len(read_names)
    name_to_idx = {name: idx for idx, name in enumerate(read_names)}

    if use_memmap:
        tmp_dir = tempfile.mkdtemp()
        matrix_file = memmap_filename or os.path.join(tmp_dir, "adj_matrix.dat")
        matrix_np = np.memmap(matrix_file, dtype='float64', mode='w+', shape=(m_size, m_size))
        matrix_np.fill(-1.0)
    else:
        matrix_np = np.full((m_size, m_size), -1.0, dtype=np.float64)

    pairs = [(read_names[i], read_names[j]) for i in range(m_size) for j in range(i + 1, m_size)]
    chunk_size = len(pairs) // (n_processes * 4) + 1
    chunks = [pairs[i:i + chunk_size] for i in range(0, len(pairs), chunk_size)]

    with Pool(processes=n_processes) as pool:
        func = partial(compute_chunk, edges=edges, data=data, snp_pos=snp_pos, R=R,
                       only_with_common_snip=only_with_common_snip)
        for chunk_result in tqdm(pool.imap_unordered(func, chunks), total=len(chunks), desc="Building matrix"):
            for a, b, d in chunk_result:
                i, j = name_to_idx[a], name_to_idx[b]
                matrix_np[i, j] = d
                matrix_np[j, i] = d

    matrix_np[0, :] = -1
    matrix_np[:, 0] = -1


    G = nx.Graph()
    total_possible_edges = (m_size * (m_size - 1)) // 2

    for idx in range(m_size):
        G.add_node(idx)

    with tqdm(total=total_possible_edges, desc="Adding edges") as pbar:
        for i in range(m_size):
            for j in range(i + 1, m_size):
                weight = matrix_np[i, j]
                if weight > weight_threshold:
                    G.add_edge(i, j, weight=weight)
                    #G.add_edge(name_to_idx[read_names[i]], name_to_idx[read_names[j]], weight=weight)
                    #G.add_edge(read_names[i], read_names[j], weight=weight)
                pbar.update(1)

    if use_memmap:
        del matrix_np

    return G
"""




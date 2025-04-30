import pysam
import networkx as nx
import numpy as np
from strainy.color_bam import color_bam2
import multiprocessing
import gfapy
import matplotlib.pyplot as plt
import matplotlib as mt
import logging
import pandas as pd
from strainy.clustering.community_detection import find_communities
from strainy.clustering import build_adj_matrix as matrix
from strainy.clustering.cluster_postprocess2 import pairs
from strainy.clustering import build_data as build_data
from strainy.params import StRainyArgs, init_global_args_storage
from strainy.params import *
import networkit as nk
from strainy.reversed import minimize
import strainy.graph_operations.gfa_ops as gfa_ops
from strainy.flye_consensus import FlyeConsensus
from strainy.overlap import StrainyOverlap2
from strainy.assembly import StrainyAssembly
from multiprocessing import Pool, cpu_count
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
import os, gzip, subprocess
from Bio import SeqIO
from Bio import SeqIO
import subprocess
import gzip
import os

logger = logging.getLogger()

def all_main(args):
    #logging.info(StRainyArgs)
    logger.info("Starting phasing")
    dirs = ("%s/vcf/" % StRainyArgs().output_intermediate,
            "%s/clusters/" % StRainyArgs().output_intermediate,
            "%s/bam/merged" % StRainyArgs().output_intermediate,
            "%s/bam/clusters" % StRainyArgs().output_intermediate,
            "%s/flye_inputs" % StRainyArgs().output_intermediate,
            "%s/graphs" % StRainyArgs().output_intermediate,
            "%s/flye_outputs/asm" % StRainyArgs().output_intermediate
)
    debug_dirs = ("%s/graphs/" % StRainyArgs().output_intermediate,
                  "%s/adj_M/" % StRainyArgs().output_intermediate
    )

    for dir in dirs:
        os.makedirs(dir, exist_ok=True)
    logger.info("### Reading SNPs...")
    snp_pos = build_data.read_snp2(StRainyArgs().snp, StRainyArgs().bam, StRainyArgs().AF)
    input_graph = gfapy.Gfa.from_file(args.gfa_ref)
    edges = input_graph.segment_names
    #edges = ["edge_192","edge_188"]
    #edges=["edge_457","edge_466","edge_467","edge_468","edge_469","edge_470","edge_471","edge_472","edge_473"]
    #print(edges)
    #unphased_edges = []
    for edge in edges:
        try:
            snp_pos[edge]
        except KeyError:
            snp_pos[edge] = []
            #unphased_edges.append(edge)
    logger.info("### Reading Reads...")
    data = build_data.read_bam2_parallel(StRainyArgs().bam, edges, snp_pos, min_mapping_quality, min_base_quality, min_al_len,
                                de_max[StRainyArgs().mode],n_processes=StRainyArgs().threads)
    cl = build_data.clusters(data)


    logger.info("### Creating connection graph...")
    g="nx" #"nk" or nx
    if g=="nx":
        m = matrix.build_adj_matrix2(edges, data, snp_pos, I, StRainyArgs().bam, 0)
        logger.info("matrix done")
        m = matrix.remove_edges(m, 0)
        m.columns = range(0, len(list(data.keys())))
        m.index = range(0, len(list(data.keys())))
        G = gfa_ops.from_pandas_adjacency_notinplace(matrix.change_w(m.transpose(), 0))
        print(G)
        """
        #TODO add NA treathment
        reads = sorted(set(cl.loc[cl["Cluster"] == "NA", "ReadName"].values))
        data_filtered = {key: data[key] for key in reads}
        m = matrix.build_adj_matrix2(edges, data_filtered, snp_pos, I, StRainyArgs().bam, 0, False)
        m = matrix.remove_edges(m, 0)
        m.columns = range(0, len(list(data_filtered.keys())))
        m.index = range(0, len(list(data_filtered.keys())))
        G = gfa_ops.from_pandas_adjacency_notinplace(matrix.change_w(m.transpose(), 0))
        indexes = cl[cl['ReadName'].isin(reads)].index.values.tolist()
        """
        logger.info("### Searching clusters...")
        cluster_membership = find_communities(G)

        clN = 0
        uncl = 0
        cl_exist = []
        for value in set(cluster_membership.values()):
            group = [k for k, v in cluster_membership.items() if v == value]
            if len(group) > 3:
                clN = clN + 1
                cl.loc[group, 'Cluster'] = value
                cl_exist.append(value)
            else:
                uncl = uncl + 1
        print("UNCLUSTERED:" + str(uncl))

      
    #TODO check clustering
    if g == "nk":
        G1=matrix.build_graph_from_data_networKit(data, snp_pos, edges, 0, only_with_common_snip=True,
                                    weight_threshold=0, n_processes=StRainyArgs().threads, chunk_size=10000)
        logger.info("### Searching clusters...")
        cluster_membership1 = nk.community.detectCommunities(G1, algo=nk.community.PLP(G1))
        clN = 0
        uncl = 0
        cl_exist = []
        for value in cluster_membership1.getSubsetIds():
            group = list(cluster_membership1.getMembers(value))
            if len(group) > 3:
                clN = clN + 1
                cl.loc[group, 'Cluster'] = value
                cl_exist.append(value)
            else:
                uncl = uncl + 1
        print("UNCLUSTERED:" + str(uncl))


    cl.to_csv(
        "%s/clusters/clusters_before_splitting_%s_%s.csv" % (StRainyArgs().output_intermediate, I, StRainyArgs().AF))
    color_bam2(edges)


    empty_consensus_dict = {}
    default_manager = multiprocessing.Manager()
    shared_flye_consensus = FlyeConsensus(StRainyArgs().bam, args.fasta_ref, StRainyArgs().threads,
                                          empty_consensus_dict, default_manager)
    reversed_edges=minimize(input_graph)["neg"]





    logger.info("### Creating overlap graph...")
    ov_nodes = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))
    ov_edges = pairs(edges, cl, shared_flye_consensus, data, snp_pos, reversed_edges, only_with_common_snip=True)
    over2 = StrainyOverlap2(ov_nodes, ov_edges)
    over2.remove_transitive()
    merged = over2.merge_unbrunching()
    for pair in merged:
        cl.loc[cl["Cluster"] == int(list(pair)[1]), "Cluster"] = int(list(pair)[0])
    over2.vis()




    logger.info("### Creating gfa graph...")
    final_graph = StrainyAssembly(over2, cl)
    final_graph.write("%s/final.gfa" % (StRainyArgs().output))
    #TODO add simplification
    color_bam2(edges)
    logger.info("Done")



if __name__ == "__main__":
    all_main()

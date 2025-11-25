import pysam
import os
import gfapy
import logging
import networkx as nx
import networkit as nk
import matplotlib as mt
from tqdm import tqdm
from lib2to3.pgen2.tokenize import group
from strainy.color_bam import color_bam2
from strainy.clustering.community_detection import find_communities
from strainy.clustering import build_adj_matrix as matrix
from strainy.clustering.cluster_postprocess import pairs
from strainy.clustering import build_data as build_data
from strainy.params import *
from strainy.flye_consensus import FlyeConsensus
from strainy.assembly import StrainyAssembly
logger = logging.getLogger()

def all_main(args):
    logger.info("Starting phasing")
    dirs = ("%s/vcf/" % StRainyArgs().output_intermediate,
            "%s/clusters/" % StRainyArgs().output_intermediate,
            "%s/bam/merged" % StRainyArgs().output_intermediate,
            "%s/bam/clusters" % StRainyArgs().output_intermediate,
            "%s/flye_inputs" % StRainyArgs().output_intermediate,
            "%s/graphs" % StRainyArgs().output_intermediate,
            "%s/flye_outputs/asm" % StRainyArgs().output_intermediate)
    debug_dirs = ("%s/graphs/" % StRainyArgs().output_intermediate,
                  "%s/adj_M/" % StRainyArgs().output_intermediate
    )
    for dir in dirs:
        os.makedirs(dir, exist_ok=True)

    logger.info("### Reading SNPs...")
    snp_pos = build_data.read_snp2(StRainyArgs().snp, StRainyArgs().bam, StRainyArgs().AF)
    input_graph = gfapy.Gfa.from_file(args.gfa_ref)
    edges = input_graph.segment_names

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


    G = matrix.build_graph_from_data_networKit(data, snp_pos, edges, 0, only_with_common_snip=True,  # TODO
                                                    weight_threshold=0, n_processes=StRainyArgs().threads,
                                                    chunk_size=10000)

    logger.info("### Creating connection graph...")
    clusters = search_clusters(G,data,cl)

    #TODO check NA cluster
    clusters = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))
    clusters_repeats = clusters

    read_names = list(data.keys())
    name_to_idx = {name: idx for idx, name in enumerate(read_names)}

    while clusters_repeats:
        clusters_repeats=iterate(clusters_repeats,cl,name_to_idx, G)

    cl.to_csv(
        "%s/clusters/clusters_before_splitting_%s_%s.csv" % (StRainyArgs().output_intermediate, I, StRainyArgs().AF))
    #cl=pd.read_csv("%s/clusters/clusters_before_splitting_%s_%s.csv" % (StRainyArgs().output_intermediate, I, StRainyArgs().AF),dtype=str, index_col=0)
    cl = cl.fillna(10000)
    сlusters = set(cl["Cluster"].values)
    color_bam2(edges)
    final_graph = StrainyAssembly(cl)
    final_graph.write("%s/final.gfa" % (StRainyArgs().output))
    logger.info("Done")


def search_clusters(G,data, cl):
    logger.info("### Searching clusters...")
    if G.numberOfEdges()==0:
        return []
    cluster_membership = nk.community.detectCommunities(G)
    clN = 0
    uncl =0
    cl_exist = sorted(set(cl.loc[cl["Cluster"] != "NA","Cluster"].values))
    new_clusters = []

    read_names = list(data.keys())
    name_to_idx = {name: idx for idx, name in enumerate(read_names)}

    for value in cluster_membership.getSubsetIds():
        group = list(cluster_membership.getMembers(value))
        names = [key for key, value in name_to_idx.items() if value in group]
        if len(group) > 3:
            new_cl_id=value
            while new_cl_id in cl_exist:
                new_cl_id = new_cl_id + 1
            clN = clN + 1
            cl.loc[cl['ReadName'].isin(names), 'Cluster'] = new_cl_id
            cl_exist.append(new_cl_id)
            new_clusters.append(new_cl_id)
        else:
            uncl = uncl + 1
    return new_clusters


def iterate(clusters_repeats,cl,name_to_idx, G):
    cl_exist = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))
    new_clusters=[]
    for cluster in clusters_repeats:
        cluster_reads = cl.loc[cl['Cluster'] == cluster]['ReadName'].values
        nodes = [name_to_idx[read] for read in cluster_reads]
        G_cl = nk.graphtools.subgraphFromNodes(G, nodes, compact=False)
        clN = 0
        uncl = 0
        cl_m = nk.community.detectCommunities(G_cl)
        if len(cl_m.getSubsetIds())>1:
            for value in cl_m.getSubsetIds():
                group = list(cl_m.getMembers(value))
                names = [key for key, value in name_to_idx.items() if value in group]
                if len(group) > 3:
                    new_cl_id = value
                    while new_cl_id in cl_exist:
                        new_cl_id = new_cl_id + 1
                    clN = clN + 1
                    cl.loc[cl['ReadName'].isin(names), 'Cluster'] = new_cl_id
                    cl_exist.append(new_cl_id)
                    new_clusters.append(new_cl_id)
                    # print("new cluster", new_cl_id)
                else:
                    uncl = uncl + 1
    return new_clusters


if __name__ == "__main__":
    all_main()


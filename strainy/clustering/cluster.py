import numpy as np
import networkx as nx
import logging
import matplotlib.pyplot as plt
import matplotlib as mt
logging.getLogger('matplotlib.font_manager').disabled = True
import multiprocessing
import pandas as pd
import pysam

from strainy.clustering.community_detection import find_communities
from strainy.clustering.cluster_postprocess import postprocess
import strainy.clustering.build_adj_matrix as matrix
import strainy.clustering.build_data as build_data
from strainy.params import *
import strainy.gfa_operations.gfa_ops as gfa_ops


logger = logging.getLogger()


def clusters_vis_stats(G, cl, clN, uncl, bam, edge, I, AF):
    cl.loc[cl['Cluster'] == 'NA', 'Cluster'] = 0
    cmap = plt.get_cmap('viridis')
    clusters=sorted(set(cl['Cluster'].astype(int)))
    cmap = cmap(np.linspace(0, 1, len(clusters)))
    colors = {}
    try:
        clusters.remove('0')
    except:
        KeyError
    colors[0] = "#505050"
    i = 0

    for cluster in clusters:
        colors[cluster] = mt.colors.to_hex(cmap[i])
        i = i + 1

    for index in cl.index:
        cl.loc[index, 'Color'] = colors[int(cl.loc[index, 'Cluster'])]
        G.remove_edges_from(list(nx.selfloop_edges(G)))
    try:
        nx.draw(G, nodelist=G.nodes(), with_labels=True, width=0.03, node_size=10, font_size=10,node_color=cl['Color'])
    except AttributeError:  #incompatability with scipy < 1.8
        pass

    ln = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4]
    cov = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[6]
    plt.suptitle(str(edge) + " coverage:" + str(cov) + " length:" + str(ln) + " clN:" + str(clN))
    plt.savefig("%s/graphs/graph_%s_%s_%s.png" % (StRainyArgs().output_intermediate, edge, I, AF), format="PNG", dpi=300)
    plt.close()

    # Calculate statistics
    logger.debug("Summary for: " + edge)
    logger.debug("Clusters found: " + str(clN))
    logger.debug("Reads unclassified: " + str(uncl))
    logger.debug("Number of reads in each cluster: ")
    logger.debug(cl['Cluster'].value_counts(dropna=False))


def cluster(i, flye_consensus):
    edge = StRainyArgs().edges_to_phase[i]
    Rcl=StRainyArgs().Rcl
    AF=StRainyArgs().AF
    R=Rcl/2
    logger.info("### Reading SNPs...")
    SNP_pos = build_data.read_snp(StRainyArgs().snp, edge, StRainyArgs().bam, AF)


    logger.info("### Reading Reads...")

    data = build_data.read_bam(StRainyArgs().bam, edge, SNP_pos, min_mapping_quality,min_base_quality, min_al_len, de_max[StRainyArgs().mode])
    cl = pd.DataFrame(columns=['ReadName', 'Cluster', 'Start'])
    for key, value in data.items():
        row = pd.DataFrame({'ReadName':[key], 'Cluster':['NA'], 'Start':[value['Start']]})
        cl = pd.concat([cl, row])
    cl = cl.reset_index(drop=True)


    total_coverage = 0
    edge_length = len(build_data.read_fasta_seq(StRainyArgs().fa, edge))
    num_reads = len(data)
    for read in data:
        total_coverage += data[read]["End"] - data[read]["Start"]
    mean_edge_cov = total_coverage // edge_length
    logger.debug(f"num reads: {num_reads}. edge length: {edge_length}, coverage: {mean_edge_cov}")

    if num_reads == 0:
        return
    if len(SNP_pos) == 0:
        #data = read_bam(StRainyArgs().bam, edge, SNP_pos, min_mapping_quality, min_al_len, de_max[StRainyArgs().mode])
        cl = pd.DataFrame(columns=['ReadName', 'Cluster', 'Start'])
        for key, value in data.items():
            row = pd.DataFrame({'ReadName':[key], 'Cluster':['NA'], 'Start':[value['Start']]})
            cl = pd.concat([cl, row])
        cl = cl.reset_index(drop=True)

        cl['Cluster'] = 1
        cl.to_csv("%s/clusters/clusters_%s_%s_%s.csv" % (StRainyArgs().output_intermediate, edge, I, AF))
        return

    #CALCULATE DISTANCE and ADJ MATRIX
    logger.info("### Calculatind distances/Building adj matrix...")
    #try:
    #    m = pd.read_csv("%s/adj_M/adj_M_%s_%s_%s.csv" % (StRainyArgs().output_intermediate, edge, I, AF), index_col='ReadName')
    #except FileNotFoundError:
    m = matrix.build_adj_matrix(cl, data, SNP_pos, I, StRainyArgs().bam, edge, R)
    if StRainyArgs().debug:
        m.to_csv("%s/adj_M/adj_M_%s_%s_%s.csv" % (StRainyArgs().output_intermediate, edge, I, AF))
    logger.info("### Removing overweighed egdes...")
    m = matrix.remove_edges(m, R)
    # BUILD graph and find clusters
    logger.info("### Creating graph...")
    m1 = m
    m1.columns = range(0,len(cl['ReadName']))
    m1.index=range(0,len(cl['ReadName']))
    G = gfa_ops.from_pandas_adjacency_notinplace(matrix.change_w(m.transpose(), R))
    logger.info("### Searching clusters...")
    cluster_membership = find_communities(G)
    clN = 0
    uncl = 0

    for value in set(cluster_membership.values()):
        group = [k for k, v in cluster_membership.items() if v == value]
        if len(group) > 3:
            clN = clN + 1
            cl.loc[group, 'Cluster'] = value
            #cl['Cluster'][group] = value
        else:
            uncl = uncl + 1

    logger.info(str(clN)+" clusters found")
    if StRainyArgs().debug:
        cl.to_csv("%s/clusters/clusters_before_splitting_%s_%s_%s.csv" % (StRainyArgs().output_intermediate, edge, I, AF))

    cl.loc[cl['Cluster'] == 'NA', 'Cluster'] = UNCLUSTERED_GROUP_N
    if clN != 0:
        logger.info("### Cluster post-processing...")
        cl = postprocess(StRainyArgs().bam, cl, SNP_pos, data, edge, R,Rcl, I, flye_consensus,mean_edge_cov)
    else:
        counts = cl['Cluster'].value_counts(dropna=False)
        cl = cl[~cl['Cluster'].isin(counts[counts < 6].index)]
    #clN = len(set(cl.loc[cl['Cluster']!='NA']['Cluster'].values))
    logger.info(str(clN) + " clusters after post-processing")
    cl.to_csv("%s/clusters/clusters_%s_%s_%s.csv" % (StRainyArgs().output_intermediate, edge, I, AF))
    
    if StRainyArgs().debug:
        logger.info("### Graph viz...")
        clusters_vis_stats(G, cl, clN,uncl, StRainyArgs().bam, edge, I, AF)

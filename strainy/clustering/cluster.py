import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mt
import multiprocessing

from strainy.clustering.community_detection import find_communities
from strainy.clustering.cluster_postprocess import postprocess
from strainy.clustering.build_adj_matrix import *
from strainy.clustering.build_data import *
from strainy.params import *


logger = logging.getLogger()


def clusters_vis_stats(G, cl, clN, uncl, SNP_pos, bam, edge, I, AF):
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

    [G.remove_node(i) for i in set(G.nodes) if i not in set(cl['ReadName'])]

    try:
        nx.draw(G, nodelist=G.nodes(), with_labels=False, width=0.03, node_size=10, font_size=5,node_color=cl['Color'])
    except:
        nx.draw(G, nodelist=G.nodes(), with_labels=False, width=0.03, node_size=10, font_size=5)

    ln = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4]
    cov = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[6]
    plt.suptitle(str(edge) + " coverage:" + str(cov) + " length:" + str(ln) + " clN:" + str(clN))
    plt.savefig("%s/graphs/graph_%s_%s_%s.png" % (StRainyArgs().output, edge, I, AF), format="PNG", dpi=300)
    plt.close()

    # Calculate statistics
    logger.debug("Summary for: " + edge)
    logger.debug("Clusters found: " + str(clN))
    logger.debug("Reads unclassified: " + str(uncl))
    logger.debug("Number of reads in each cluster: ")
    logger.debug(cl['Cluster'].value_counts(dropna=False))



def cluster(i, flye_consensus):
    edge = StRainyArgs().edges[i]
    logger.info("### Reading SNPs...")
    SNP_pos = read_snp(StRainyArgs().snp, edge, StRainyArgs().bam, AF)

    logger.info("### Reading Reads...")
    data = read_bam(StRainyArgs().bam, edge, SNP_pos, clipp, min_mapping_quality, min_al_len, de_max[StRainyArgs().mode])
    cl=pd.DataFrame(data={'ReadName': data.keys()})
    logger.debug(str(len(cl['ReadName'])) + " reads found")
    cl['Cluster'] = 'NA'
    if len(cl['ReadName'])==0:
        return
    if len(SNP_pos)==0:
        #data = read_bam(StRainyArgs().bam, edge, SNP_pos, clipp, min_mapping_quality, min_al_len, de_max[StRainyArgs().mode])
        cl = pd.DataFrame(data={'ReadName': data.keys()})
        cl['Cluster'] = 1
        cl.to_csv("%s/clusters/clusters_%s_%s_%s.csv" % (StRainyArgs().output, edge, I, AF))
        return

    #CALCULATE DISTANCE and ADJ MATRIX
    logger.info("### Calculatind distances/Building adj matrix...")
    #try:
    #    m = pd.read_csv("%s/adj_M/adj_M_%s_%s_%s.csv" % (StRainyArgs().output, edge, I, AF), index_col='ReadName')
    #except FileNotFoundError:
    m = build_adj_matrix(cl, data, SNP_pos, I, StRainyArgs().bam, edge, R)
    m.to_csv("%s/adj_M/adj_M_%s_%s_%s.csv" % (StRainyArgs().output, edge, I, AF))

    logger.info("### Removing overweighed egdes...")
    m = remove_edges(m, R)

    # BUILD graph and find clusters
    logger.info("### Creating graph...")
    m1 = m
    m1.columns = range(0,len(cl['ReadName']))
    m1.index=range(0,len(cl['ReadName']))
    G = nx.from_pandas_adjacency(change_w(m.transpose(), R))
    logger.info("### Searching clusters...")
    cluster_membership = find_communities(G)
    clN = 0
    uncl = 0

    for value in set(cluster_membership.values()):
        group = [k for k, v in cluster_membership.items() if v == value]
        if len(group) > 3:
            clN = clN + 1
            cl['Cluster'][group] = value
        else:
            uncl = uncl + 1

    logger.info(str(clN)+" clusters found")
    cl.to_csv("%s/clusters/clusters_before_splitting_%s_%s_%s.csv" % (StRainyArgs().output, edge, I, AF))

    cl.loc[cl['Cluster'] == 'NA', 'Cluster'] = unclustered_group_N
    if clN != 0:
        logger.info("### Cluster post-processing...")
        cl = postprocess(StRainyArgs().bam, cl, SNP_pos, data, edge, R, I, flye_consensus)
    else:
        counts = cl['Cluster'].value_counts(dropna=False)
        cl = cl[~cl['Cluster'].isin(counts[counts < 6].index)]
    clN = len(set(cl.loc[cl['Cluster']!='NA']['Cluster'].values))
    logger.info(str(clN) + " clusters after post-processing")
    cl.to_csv("%s/clusters/clusters_%s_%s_%s.csv" % (StRainyArgs().output, edge, I, AF))
    logger.info("### Graph viz...")

    clusters_vis_stats (G,cl, clN,uncl,SNP_pos, StRainyArgs().bam, edge, I, AF)
    cl.to_csv("%s/clusters/clusters_%s_%s_%s.csv" % (StRainyArgs().output, edge, I, AF))

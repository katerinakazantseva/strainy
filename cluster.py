import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mt

from community_detection import find_communities
from cluster_postprocess import postprocess
from build_adj_matrix import *
from build_data import *
from params import *


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
        cl.loc[index, 'Color'] = colors[cl.loc[index, 'Cluster']]
    G.remove_edges_from(list(nx.selfloop_edges(G)))

    nx.draw(G, nodelist=G.nodes(), with_labels=False, width=0.03, node_size=10, font_size=5,
            node_color=cl['Color'])

    ln = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4]
    cov = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[6]
    plt.suptitle(str(edge) + " coverage:" + str(cov) + " length:" + str(ln) + " clN:" + str(clN))
    plt.savefig("output/graphs/graph_%s_%s_%s.png" % (edge, I, AF), format="PNG", dpi=300)
    plt.close()

    # Calculate statistics
    print("Summary for: " + edge)
    print("Clusters found: " + str(clN))
    print("Reads unclassified: " + str(uncl))
    print("Number of reads in each cluster: ")
    print(cl['Cluster'].value_counts(dropna=False))

    stats = open('output/stats.txt', 'a')
    stats.write(edge + "\t" + str(ln) + "\t" + str(cov) + "\t" + str(len(cl['ReadName'])) + "\t" + str(
        len(SNP_pos)) + "\t" + str(clN) + "\t" + str(uncl) + "\n")
    stats.close()


def cluster(params):
    # params = #i, consensus_dict)
    i, flye_consensus = params
    edge=edges[i]
    print("### Reading SNPs...")
    SNP_pos = read_snp(snp, edge, bam, AF)
    if len(SNP_pos) == 0:
        return

    # READ READS AND POSITIONS
    print("### Reading Reads...")
    data = read_bam(bam,edge,SNP_pos,clipp,min_mapping_quality,min_al_len,de_max)
    cl = pd.DataFrame(data={'ReadName': data.keys()})
    cl['Cluster'] = 'NA'
    print(str(len(cl['ReadName']))+" reads found")
    if len(cl['ReadName']) == 0:
        return

    # CALCULATE DISTANCE and ADJ MATRIX
    print("### Calculatind distances/Building adj matrix...")
    m = build_adj_matrix(cl, data, SNP_pos, I, bam, edge, R)
    #m = pd.read_csv("output/adj_M/adj_M_%s_%s_%s.csv" % (edge, I, AF), index_col='ReadName')
    m.to_csv("output/adj_M/adj_M_%s_%s_%s.csv" % (edge, I, AF))
    print("### Removing overweighed egdes...")
    m = remove_edges(m, R)

    # BUILD graph and find clusters
    print("### Creating graph...")
    m1 = m
    m1.columns = range(0,len(cl['ReadName']))
    m1.index=range(0,len(cl['ReadName']))
    G = nx.from_pandas_adjacency(change_w(m.transpose(), R))
    print("### Searching clusters...")
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

    print(str(clN)+" clusters found")
    cl.to_csv("output/clusters/clusters_before_splitting_%s_%s_%s.csv" % (edge, I, AF))
    print("### Cluster post-processing...")
    cl.loc[cl['Cluster'] == 'NA', 'Cluster'] = 1000000
    if clN != 0:
        cl = postprocess(bam, cl, SNP_pos, data, edge, R, I, flye_consensus)
    clN = len(set(cl.loc[cl['Cluster']!='NA']['Cluster'].values))
    print(str(clN) + " clusters after post-processing")
    cl.to_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF))
    print("### Graph viz...")
    clusters_vis_stats(G, cl, clN, uncl, SNP_pos, bam, edge, I, AF)


stats = open('output/stats.txt', 'a')
stats.write("Edge" + "\t" + "Len" + "\t" + "Coverage" + "\t" + "ReadsN" + "\t"+"SNPN"+"\t"+"ClustersN"+"\t"+"UnclusteredRN"+"\n")
stats.close()
























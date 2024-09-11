import logging
from strainy.graph_operations import gfa_ops
from strainy.unitig_statistics import utg_stats
from strainy.clustering import build_data
from strainy.params import *

"""
This contains functions for operation with assembly graph:
1. add_child_edge: adds a child unitig with the same sequence as the father unitig or with the given sequence 
(aff full path and othr clusters"
2. add_path_links: Adds link ("full path)
3.change_cov: recalculate coverage
4.change_sec" recalculate sequence
"""




logger = logging.getLogger()



def add_child_edge(edge, clN, g, cl, left, right, cons, flye_consensus, change_seq=True, insertmain=True):
    """
    Adds a child unitig with the same sequence as the parental unitig or with the given sequence
    """
    ##TODO if cons provided change_seq=True (provide seq not consensus)
    ##TODO make separare function to add gfa edge and move to gfa_ops
    consensus = flye_consensus.flye_consensus(clN, edge, cl)
    consensus_start = consensus["start"]
    cons_length_diff = len(consensus["consensus"]) - (consensus["end"] - consensus["start"])
    logger.debug(f'Consensus length difference: {cons_length_diff}')
    if consensus_start > left and insertmain == True:
        main_seq = g.try_get_segment(edge)
        insert = main_seq.sequence[left:consensus_start]
        seq = str(consensus["consensus"])[0: right - consensus_start + cons_length_diff + 1]
        seq = insert + seq
    else:
        seq = str(consensus["consensus"])[left - consensus_start: right - consensus_start + cons_length_diff + 1]

    new_line = gfa_ops.add_edge(g, edge, clN, round(cons[clN]["Cov"]))

    if change_seq == True:  ##TODO: move to gfa_ops.add_edge
        if len(seq) == 0:
            new_line.sequence = "A"
        else:
            new_line.sequence = seq
    else:
        new_line.sequence = g.try_get_segment("%s" % edge).sequence

    logger.debug("Unitig created  %s_%s" % (edge, clN))
    utg_stats.store_phased_unitig_info(new_line,
                                       edge,
                                       len(cons[clN]) - 7,
                                       left,
                                       right
                                       )


def add_path_links(graph, edge, paths):
    """
     Add gfa links between newly created unitigs forming "full path"
    """
    for path in paths:
        for i in range(0, len(path) - 1):
            gfa_ops.add_link(graph, f"{edge}_{path[i]}", "+", f"{edge}_{path[i + 1]}", "+", 1)



def change_cov(g, edge, cons, ln, clusters, othercl, remove_clusters):
    cov = 0
    len_cl = []
    for i in othercl:
        cov += cons[i]["Cov"] * (cons[i]["End"] - cons[i]["Start"])
        for i in range(cons[i]["Start"],cons[i]["End"]):
            len_cl.append(i)
    if (len(set(len_cl)) / ln) < parental_min_len and len(clusters)- len(othercl) != 0:
        remove_clusters.add(edge)
    cov = cov / ln
    i = g.try_get_segment(edge)
    i.dp = round(cov)
    return cov


def change_sec(g, edge, othercl, cl,SNP_pos, data, cut = True):
    temp = {}
    other_cl = cl
    for cluster in othercl:
        other_cl.loc[cl["Cluster"] == cluster, "Cluster"] = "OTHER_%s" %edge

    reference_seq = build_data.read_fasta_seq(StRainyArgs().fa, edge)
    cl_consensuns = build_data.cluster_consensuns(other_cl, "OTHER_%s" %edge, SNP_pos, data, temp, edge, reference_seq)
    i = g.try_get_segment(edge)
    seq = i.sequence
    seq = list(seq)
    for key, val in cl_consensuns["OTHER_%s" %edge].items():
        try:
            seq[int(key) - 1] = val
        except (ValueError):
            continue
    i.sequence=''.join(seq)


def strong_tail(cluster, cl, ln, data):
    count_start = None
    count_stop = None
    res = [False,False]
    reads = list(cl.loc[cl["Cluster"] == cluster, "ReadName"])
    for read in reads:
        if data[read]["Start"] < start_end_gap:
            if count_start == None:
                count_start = 0
            count_start = count_start+1
        if data[read]["End"] > ln - start_end_gap:
            if count_stop == None:
                count_stop = 0
            count_stop = count_stop + 1
    if count_start!= None and count_start > strong_cluster_min_reads :
        res[0] = True
    if  count_stop!= None and count_stop > strong_cluster_min_reads:
        res[1] = True
    return (res)


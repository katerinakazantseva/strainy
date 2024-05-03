import networkx as nx
import pygraphviz as gv
import re
import gfapy
from collections import Counter, deque, defaultdict
import pandas as pd
import pickle
import logging
import multiprocessing
import shutil
import pysam
import time
import traceback
import csv

import strainy.clustering.build_adj_matrix as matrix
import strainy.clustering.cluster_postprocess as postprocess
import strainy.simplification.simplify_links as smpl
import strainy.gfa_operations.gfa_ops as gfa_ops
from strainy.flye_consensus import FlyeConsensus
import strainy.clustering.build_data as build_data
from strainy.params import *
from strainy.logging import set_thread_logging
from strainy.reports.strainy_stats import strain_stats_report
from strainy.reports.call_variants import produce_strainy_vcf
from strainy.preprocessing import gfa_to_fasta

logger = logging.getLogger()

def format_rounding(number):
    n = abs(number)
    if n == 0:
        return '0.000'
    if n < 1:
        # Find the first non-zero digit.
        # We want 3 digits, starting at that location.
        s = f'{n:.99f}'
        index = re.search('[1-9]', s).start()
        return s[:index + 3]
    else:
        # We want 2 digits after decimal point.
        return str(round(n, 2))
    

def write_phased_unitig_csv():
    columns = [
        'Strain_unitig',
        'Reference_unitig',
        'Length',
        'Coverage',
        'Abundance_Ratio',
        '#SNP',
        'SNP_density',
        'Start_positioin',
        'End_position'
        ]
    
    with open(StRainyArgs().phased_unitig_info_table_path, 'w') as f:
        write = csv.writer(f, delimiter='\t')

        write.writerow(columns)
        write.writerows(list(StRainyArgs().phased_unitig_info_table.values()))


def write_reference_unitig_csv():
    columns=[
        'Reference_unitig',
        'Length',
        'Coverage',
        'SNP_density',
        'Is_processed',
        'Is_phased'
        ]
    
    with open(StRainyArgs().reference_unitig_info_table_path, 'w') as f:
        write = csv.writer(f, delimiter='\t')
        
        write.writerow(columns)
        write.writerows(list(StRainyArgs().reference_unitig_info_table.values()))


def store_phased_unitig_info(strain_unitig, reference_unitig, n_SNPs, start, end):
    reference_coverage = round(float(pysam.samtools.coverage("-r",
                                                             reference_unitig,
                                                             StRainyArgs().bam,
                                                             "--no-header").
                                                             split()[6]))
    # # Log the information to std output
    # logger.info(f'== == Inserted Strain unitig: {strain_unitig.name} == == ')
    # logger.info(f'\t\t Reference unitig: {reference_unitig}')
    # logger.info(f'\t\t Length: {strain_unitig.length} bp')
    # logger.info(f'\t\t Coverage: {strain_unitig.dp}')
    # logger.info(f'\t\t Abundance Ratio: {round(100 * strain_unitig.dp // reference_coverage)}%')
    # logger.info(f'\t\t #SNP: {n_SNPs}')
    # logger.info(f'\t\t SNP density: {format_rounding(n_SNPs / strain_unitig.length)}')
    # logger.info(f'\t\t Start position: {start}')
    # logger.info(f'\t\t End position: {end}\n\n')

    try:
        abundance_ratio =  round(100 * strain_unitig.dp // reference_coverage)
    except ZeroDivisionError:
        abundance_ratio = 0

    StRainyArgs().phased_unitig_info_table[strain_unitig.name] = [
        strain_unitig.name,
        reference_unitig,
        strain_unitig.length,
        strain_unitig.dp,
        abundance_ratio,
        n_SNPs,
        format_rounding(n_SNPs / strain_unitig.length),
        start,
        end
        ]
    

def store_reference_unitig_info(ref_coverage):
    graph = gfapy.Gfa.from_file(StRainyArgs().gfa)
    phased_unitig_df = pd.read_csv(StRainyArgs().phased_unitig_info_table_path, sep='\t')
    counter = Counter(list(phased_unitig_df['Reference_unitig']))
    for reference_unitig in graph.segments:

        # Number of phased unitigs created from this reference unitig
        n_phased_unitigs = counter[reference_unitig.name]
        # Number of SNPs
        n_SNPs = len(
            build_data.read_snp(
                StRainyArgs().snp,
                reference_unitig.name,
                StRainyArgs().bam,
                StRainyArgs().AF
                )
            )
        StRainyArgs().reference_unitig_info_table[reference_unitig.name] = [
            reference_unitig.name,
            reference_unitig.length,
            ref_coverage[reference_unitig.name],
            format_rounding(n_SNPs / reference_unitig.length),
            reference_unitig.name in StRainyArgs().edges_to_phase,
            n_phased_unitigs > 1
        ]


def add_child_edge(edge, clN, g, cl, left, right, cons, flye_consensus, change_seq=True, insertmain=True):
    """
    The function creates unitigs in the gfa graph
    """
    ##TODO make separare function to add gfa edge and move to gfa_ops
    consensus = flye_consensus.flye_consensus(clN, edge, cl)
    consensus_start = consensus["start"]
    cons_length_diff = len(consensus["consensus"]) - (consensus["end"] - consensus["start"])
    logger.debug(f'Consensus length difference: {cons_length_diff}')
    if consensus_start > left and insertmain==True:
        main_seq = g.try_get_segment(edge)
        insert = main_seq.sequence[left:consensus_start]
        seq = str(consensus["consensus"])[0 : right - consensus_start + cons_length_diff + 1]
        seq = insert+seq
    else:
        seq = str(consensus["consensus"])[left - consensus_start : right - consensus_start + cons_length_diff + 1]

    g.add_line("S\t%s_%s\t*" % (edge, clN))
    new_line = g.try_get_segment("%s_%s" % (edge, clN))
    new_line.name = str(edge) + "_" + str(clN)
    new_line.sid = str(edge) + "_" + str(clN)
    new_line.dp = round(cons[clN]["Cov"])  # TODO: what to do with coverage?
    #remove_zeroes.append("S\t%s_%s\t*" % (edge, clN))
    if change_seq==True:
        if len(seq) == 0:
            new_line.sequence = "A"
        else:
            new_line.sequence = seq
    else:
        new_line.sequence = g.try_get_segment("%s" % edge).sequence

    logger.debug("Unitig created  %s_%s" % (edge, clN))

    store_phased_unitig_info(new_line,
                    edge,
                    len(cons[clN]) - 7,
                    left,
                    right
                    )
    

def build_paths_graph(cons, full_paths_roots, full_paths_leafs, cluster_distances):
    """
    Create an "overlap" graph for clusters within a unitig, based on flye distance
    """
    M = cluster_distances
    G = gfa_ops.from_pandas_adjacency_notinplace(M, create_using = nx.DiGraph)
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    try:
        G.remove_node(0)
    except:
        pass
    path_remove = []
    node_remove = []
    for node in full_paths_leafs:
        neighbors = list(full_paths_leafs)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G, node, neighbor, cutoff = 2):
                if len(n_path) == 2:
                    node_remove.append(neighbor)

    for node in full_paths_roots:
        neighbors = list(full_paths_roots)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G,  neighbor,node, cutoff = 2):
                if len(n_path) == 2:
                    node_remove.append(neighbor)
    G = remove_nested(G, cons)
    for node in node_remove:
        try:
            G.remove_node(node)
            logger.debug("REMOVE " + str(node))
            full_paths_roots.remove(node)
            full_paths_leafs.remove(node)
        except:
            continue

    for node in G.nodes():
        neighbors = nx.all_neighbors(G, node)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G, node, neighbor, cutoff = 3):
                if len(n_path) == 3:
                    path_remove.append(n_path)

    for n_path in path_remove:
        try:
            G.remove_edge(n_path[0], n_path[1])
        except:
            continue
    return (G)


def remove_nested(G, cons):
    """
     Disconnect "nested" clusters from the parent cluster
    """
    nodes = list(G.nodes())
    for node in nodes:
        try:
            neighbors = nx.all_neighbors(G, node)
            for neighbor in list(neighbors):
                if cons[node]["Start"] < cons[neighbor]["Start"] and cons[node]["End"] > cons[neighbor]["End"]:
                    try:
                        G.remove_edge(node, neighbor)
                        G.remove_edge(neighbor,node)
                        logger.debug("REMOVE NESTED" + str(neighbor))

                    except:
                        continue
        except:
            continue
    return (G)


def paths_graph_add_vis(edge, cons, cl, full_paths_roots,
                        full_paths_leafs, full_clusters, cluster_distances):
    """
     Graph visualization function
    """
    G_vis = gfa_ops.from_pandas_adjacency_notinplace(cluster_distances,
                                                     create_using = nx.DiGraph)
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))

    try:
        G_vis.remove_node(0)
    except:
        pass

    cluster_colors = {}
    for i, row in cl.iterrows():
        if row["Cluster"] not in cluster_colors:
            cluster_colors[row["Cluster"]] = row["Color"]

    for e in G_vis.edges():
        first_cl, second_cl = e
        intersect = min(cons[first_cl]["End"], cons[second_cl]["End"]) - \
                    max(cons[first_cl]["Start"], cons[second_cl]["Start"])
        G_vis[first_cl][second_cl]["label"] = f"Ovlp:{intersect}"

    for n in G_vis.nodes():
        clust_len = cons[n]["End"] - cons[n]["Start"]
        color = cluster_colors[n]
        G_vis.nodes[n]["label"] = f"{color} len:{clust_len}"

    G_vis.add_node("Src",style = "filled",fillcolor = "gray",shape = "square")
    G_vis.add_node("Sink",style = "filled",fillcolor = "gray",shape = "square")
    for i in full_paths_roots:
        G_vis.add_edge("Src", i)

    for i in full_paths_leafs:
        G_vis.add_edge(i, "Sink")

    for i in full_clusters:
        G_vis.add_edge("Src", i)
        G_vis.add_edge(i, "Sink")

    graph_str = str(nx.nx_agraph.to_agraph(G_vis))
    graph_vis = gv.AGraph(graph_str)
    graph_vis.layout(prog = "dot") # TODO: this line may cause an error
    graph_vis.draw("%s/graphs/connection_graph_%s.png" % (StRainyArgs().output_intermediate, edge))


def find_full_paths(G, paths_roots, paths_leafs):
    paths = []
    for root in paths_roots:
        try:
            #TODO: we need to increas cutoff for longer unitigs with more clusters.
            #But this will result in the exponential number of paths. Instead, we should be
            #looking at all nodes that are reachable from both source and sink, which is linear
            paths_nx = nx.algorithms.all_simple_paths(G, root, paths_leafs, cutoff = 10)
        except:
            pass
        for path in list(paths_nx):
            paths.append(path)

    return (paths)


def add_path_links(graph, edge, paths, G):
    """
     Add gfa links between newly created unitigs forming "full path"
    """
    for path in paths:
        for i in range(0, len(path) - 1):
            gfa_ops.add_link(graph, f"{edge}_{path[i]}", "+", f"{edge}_{path[i + 1]}", "+", 1)


def add_path_edges(edge, g, cl, ln, full_paths, G, paths_roots, paths_leafs, full_clusters, cons, flye_consensus):
    """
    Add gfa nodes (unitigs) forming "full path", calculating cluster boundaries
    """
    path_cl = []
    logger.debug("Add path")
    for node in full_clusters:
        try:
            paths_roots.remove(node)
            paths_leafs.remove(node)
        except:
            pass

    for path in full_paths.copy():
        for member in path:
            if member in full_clusters:
                try:
                    full_paths.remove(path)
                except (ValueError):
                    pass
            if member in paths_leafs and path.index(member)!= len(path)-1:
                try:
                    full_paths.remove(path)
                except (ValueError):
                    pass

    for path in full_paths:
        for member in path:
            path_cl.append(member)
    cut_l_unsorted = {}
    cut_r_unsorted = {}
    for path_cluster in set(path_cl):
        cut_l_unsorted[path_cluster] = None
        cut_r_unsorted[path_cluster] = None
        if path_cluster in paths_roots and cons[path_cluster]["Start"] < start_end_gap :
            cut_l_unsorted[path_cluster] = cons[path_cluster]["Start"]
        if path_cluster in paths_leafs:
            cut_r_unsorted[path_cluster] = ln - 1
    stop_pos = {}
    for i in cut_r_unsorted.keys():
        stop_pos[i] = cons[i]["End"]
    order_by_stop_pos = list(dict(sorted(stop_pos.items(), key = lambda item: item[1])).keys())

    cut_l = {}
    cut_r = {}
    for i in order_by_stop_pos:
        cut_l[i] = cut_l_unsorted[i]
        cut_r[i] = cut_r_unsorted[i]
    Members=list(cut_l.keys())
    while Members:
        member=Members.pop(0)
        if cut_l[member] == None and cut_r[member] == None: #if the cluster does not already have boundaries, try the next one first
                member_to_q=member
                member=Members.pop(0)
                Members.insert(0,member_to_q)
        if cut_l[member] != None and (cut_r[member] == None or member in paths_leafs):
            Q = deque()
            L = []
            R = []
            for path in full_paths:
                try:
                    L.append(path[path.index(member) + 1])
                    Q.append(path[path.index(member) + 1])
                except (ValueError, IndexError):
                    continue
            visited = []
            Q = list(set(Q))
            while Q:
                n = Q.pop()
                visited.append(n)
                if n in L:
                    for path in full_paths:
                        try:
                            if path.index(n) > 0:
                                if path[path.index(n) - 1] not in visited:
                                    R.append(path[path.index(n) - 1])
                                    if path[path.index(n) - 1] not in Q:
                                        Q.append(path[path.index(n) - 1])
                        except (ValueError, IndexError):
                            continue
                else:
                    for path in full_paths:
                        try:
                            if path[path.index(n) + 1] not in visited:
                                L.append(path[path.index(n) + 1])
                                if path[path.index(n) + 1] not in Q:
                                    Q.append(path[path.index(n) + 1])

                        except (ValueError, IndexError):
                                continue
            l_borders = []
            r_borders = []
            for i in L:
                l_borders.append(int(cons[i]["Start"]))

            for i in R:
                r_borders.append(int(cons[i]["End"]))
            if member in paths_leafs:
                border = cut_r[member]
            else:
                border = max(l_borders) + (min(r_borders) - max(l_borders)) // 2
            for i in L:
                cut_l[i] = border
            for i in R:
                cut_r[i] = border
        elif cut_r[member] != None:
            for path in full_paths:
                try:
                    cut_l[path[path.index(member)+1]] = cut_r[member]
                except:
                    pass

    if None in cut_l.values():
        for member in cut_l.keys():
            if cut_l[member] == None:
                for path in full_paths:
                    try:
                        cut_l[member] = cut_r[path[path.index(member)-1]]
                    except:
                        pass
    for path_cluster in set(path_cl):
        if cut_l[path_cluster]!= cut_r[path_cluster]:
            add_child_edge(edge, path_cluster, g,  cl, cut_l[path_cluster], cut_r[path_cluster], cons, flye_consensus)
        else:
            for i in range(0,len(full_paths)):
                if path_cluster in full_paths[i]:
                    upd_path = full_paths[i]
                    upd_path.remove(path_cluster)
                    full_paths[i] = upd_path
            G.remove_node(path_cluster)

    return(path_cl)


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


def gcu_worker(edge, flye_consensus, args):
    init_global_args_storage(args)

    bam_cache = {}
    link_clusters = defaultdict(list)
    link_clusters_src = defaultdict(list)
    link_clusters_sink = defaultdict(list)
    graph_ops = []
    remove_clusters = set()

    set_thread_logging(StRainyArgs().log_transform, "gcu", multiprocessing.current_process().pid)
    logger.info("\n\n\t == == Processing unitig " + edge + " == == ")

    try:
        graph_create_unitigs(edge,
                            flye_consensus,
                            bam_cache,
                            link_clusters,
                            link_clusters_src,
                            link_clusters_sink,
                            remove_clusters,
                            graph_ops)
    except Exception as e:
        logger.error("Worker thread exception! " + str(e) + "\n" + traceback.format_exc())
        raise e
    return bam_cache, link_clusters, link_clusters_src, link_clusters_sink, graph_ops, remove_clusters


def parallelize_gcu(pool, graph_edges, flye_consensus, graph, args):
    if StRainyArgs().threads == 1:
        result_values = []
        for edge in graph_edges:
            result_values.append(gcu_worker(edge, flye_consensus, args))

    else:
        init_args = [(edge, flye_consensus, args) for edge in graph_edges]
        results = pool.starmap_async(gcu_worker, init_args, chunksize=1)
        while not results.ready():
            time.sleep(0.01)
            if not results._success:
                pool.terminate()
                raise Exception("Error in worker thread, exiting")
        result_values = results._value
        pool.close()
        pool.join()

    bam_cache = {}
    link_clusters = defaultdict(list)
    link_clusters_src = defaultdict(list)
    link_clusters_sink = defaultdict(list)
    graph_ops = []
    remove_clusters = []
    
    outputs = [bam_cache, link_clusters, link_clusters_src, link_clusters_sink, graph_ops, remove_clusters]
    # join the results of multiple threads
    for r in result_values:
        for i in range(len(r)):
            if i == len(r) - 1 or i == len(r) - 2:
                for k in r[i]:
                    outputs[i].append(k)
            else:
                for k, v in r[i].items():
                    outputs[i][k] = v
    
    # operations on the graph performed after the parallel graph_create_unitig
    # this is due to not being able to pass the graph object to threads
    for op in graph_ops:
        if op[0] == 'add_child_edge':
            add_child_edge(op[1], op[2], graph, op[3], op[4], op[5], op[6], flye_consensus, op[7], op[8])
        elif op[0] == 'add_path_edges':
            add_path_edges(op[1], graph, op[2], op[5], op[6], op[7], op[8], op[9], op[10], op[11], flye_consensus)
        elif op[0] == 'add_path_links':
            add_path_links(graph, op[1], op[2], op[3])
        
    return bam_cache, link_clusters, link_clusters_src, link_clusters_sink, set(remove_clusters), graph


def graph_create_unitigs(edge, flye_consensus, bam_cache, link_clusters,
                         link_clusters_src, link_clusters_sink, remove_clusters, graph_ops):
    """
    First part of the transformation: creation of all new unitigs from clusters obtained during the phasing stage
    """
    full_paths_roots = []
    full_paths_leafs = []
    full_paths = []
    full_clusters = []

    graph = gfapy.Gfa.from_file(StRainyArgs().gfa)

    cl = None
    try:
        cl = pd.read_csv("%s/clusters/clusters_%s_%s_%s.csv" % (StRainyArgs().output_intermediate, edge, I, StRainyArgs().AF), keep_default_na = False)
    except(FileNotFoundError, IndexError):
        logger.debug("%s: No clusters" % edge)
        clusters = []

    if cl is not None:
        SNP_pos = build_data.read_snp(StRainyArgs().snp, edge, StRainyArgs().bam, StRainyArgs().AF)
        data = build_data.read_bam(StRainyArgs().bam, edge, SNP_pos, min_mapping_quality,min_base_quality,min_al_len, de_max[StRainyArgs().mode])
        bam_cache[edge] = data

        ln = int(pysam.samtools.coverage("-r", edge, StRainyArgs().bam, "--no-header").split()[4])
        if len(cl.loc[cl["Cluster"] == 0,"Cluster"].values) > 10:
            cl.loc[cl["Cluster"] == 0, "Cluster"] = 1000000
        clusters = sorted(set(cl.loc[cl["Cluster"] != "NA","Cluster"].values))

        try:
            clusters.remove(0)
        except:
            pass

        reference_seq = build_data.read_fasta_seq(StRainyArgs().fa, edge)
        cons = build_data.build_data_cons(cl, SNP_pos, data, edge, reference_seq)

        if len(clusters) == 1:
            for cluster in clusters:
                clStart = cons[cluster]["Start"]
                clStop = cons[cluster]["End"]
                if clStart < start_end_gap and clStop > ln - start_end_gap:
                    full_paths_roots.append(cluster)
                    full_paths_leafs.append(cluster)
                consensus = flye_consensus.flye_consensus(cluster, edge, cl)
                # add_child_edge(edge, cluster, graph, cl, consensus["start"], consensus["end"], cons, flye_consensus,change_seq = False)
                graph_ops.append(['add_child_edge', edge, cluster, cl, consensus["start"], consensus["end"], cons, False, True])
            link_clusters[edge] = list(clusters)
            link_clusters_sink[edge] = list(clusters)
            link_clusters_src[edge] = list(clusters)
            remove_clusters.add(edge)

        if len(clusters) > 1:
            for cluster in clusters:
                clStart = cons[cluster]["Start"]
                clStop = cons[cluster]["End"]
                if clStart < start_end_gap and clStop > ln - start_end_gap:
                    if strong_tail(cluster, cl, ln, data)[0] == True and strong_tail(cluster, cl, ln,data)[1] == True:
                        consensus = flye_consensus.flye_consensus(cluster, edge, cl)
                        # add_child_edge(edge, cluster, graph, cl,consensus["start"], consensus["end"], cons, flye_consensus)
                        graph_ops.append(['add_child_edge', edge, cluster, cl, consensus["start"], consensus["end"], cons, True, True])
                        full_clusters.append(cluster)

                    elif strong_tail(cluster, cl, ln, data)[0] != True:
                        cons[cluster]["Start"] = cons[cluster]["Start"] + start_end_gap+1
                    else:
                        cons[cluster]["End"] = cons[cluster]["End"] - start_end_gap-1
                if clStart < start_end_gap and strong_tail(cluster, cl, ln, data)[0] == True :
                    full_paths_roots.append(cluster)
                if clStop > ln - start_end_gap and strong_tail(cluster, cl, ln, data)[1] == True:
                    full_paths_leafs.append(cluster)

            cluster_distances = postprocess.build_adj_matrix_clusters(edge, cons, cl, flye_consensus, False)
            cluster_distances = matrix.change_w(cluster_distances,0)

            G = build_paths_graph(cons, full_paths_roots, full_paths_leafs, cluster_distances.copy())

            #full_cl[edge] = full_clusters
            if StRainyArgs().debug:
                paths_graph_add_vis(edge, 
                                    cons,
                                    cl,
                                    full_paths_roots,
                                    full_paths_leafs,
                                    full_clusters,
                                    cluster_distances.copy())

            try:
                full_paths = find_full_paths(G,full_paths_roots, full_paths_leafs)
            except(ValueError):
                pass

            # add_path_edges(edge, graph, cl, data, SNP_pos, ln, full_paths, G,full_paths_roots,
                        #    full_paths_leafs,full_clusters,cons, flye_consensus)
            graph_ops.append(['add_path_edges', edge, cl, data, SNP_pos, ln, full_paths, G,full_paths_roots,
                           full_paths_leafs,full_clusters,cons])
            # add_path_links(graph, edge, full_paths, G)
            graph_ops.append(['add_path_links', edge, full_paths, G])

            othercl = list(set(clusters) - set(full_clusters) - set([j for i in full_paths for j in i]))
            if len(othercl) > 0:
                G = gfa_ops.from_pandas_adjacency_notinplace(cluster_distances.copy(), create_using = nx.DiGraph)

            close_to_full = []
            othercl_len=[cons[i]['End']-cons[i]['Start'] for i in othercl]
            othercl_sorted=[i[1] for i  in sorted(zip(othercl_len, othercl), reverse=True)]
            removed=[]
            for cluster in othercl_sorted:
                neighbors = nx.all_neighbors(G, cluster)
                A = set(neighbors)
                B = set([j for i in full_paths for j in i])
                if len(A.intersection(set(full_clusters))) > 0 or len(A.intersection(B)) > 0: #remove close-to full to avoid duplication
                    try:
                        othercl.remove(cluster)
                        close_to_full.append(cluster)
                        removed.append(cluster)
                    except (ValueError):
                        pass
                if len(A)>0 and cluster not in removed: #leave longest and remove their neighbors
                    for i in A:
                        try:
                            othercl.remove(i)
                            removed.append(i)
                            logger.debug("REMOVE " + str(cluster))
                        except (ValueError):
                            pass


            new_cov = change_cov(graph, edge, cons, ln, clusters, othercl, remove_clusters)
            if  new_cov < parental_min_coverage and len(clusters) - len(othercl) != 0 and (len(set(full_clusters))>0 or len(full_paths)>0):
                remove_clusters.add(edge)
            else:
                for cluster in othercl:
                    consensus = flye_consensus.flye_consensus(cluster, edge, cl)
                    # add_child_edge(edge, cluster, graph, cl, cons[cluster]["Start"], cons[cluster]["End"], cons, flye_consensus,insertmain=False)
                    graph_ops.append(['add_child_edge', edge, cluster, cl, cons[cluster]["Start"], cons[cluster]["End"], cons, True, False])
                remove_clusters.add(edge)

            link_clusters[edge] = list(full_clusters) + list(
                set(full_paths_roots).intersection(set([j for i in full_paths for j in i]))) + list(
                set(full_paths_leafs).intersection(set([j for i in full_paths for j in i])))
            link_clusters_src[edge] = list(full_clusters) + list(
                set(full_paths_roots).intersection(set([j for i in full_paths for j in i])))
            link_clusters_sink[edge] = list(full_clusters) + list(
                set(full_paths_leafs).intersection(set([j for i in full_paths for j in i])))

    stats = open("%s/stats_clusters.txt" % StRainyArgs().output_intermediate, "a")
    fcN = 0
    fpN = 0

    try:
        fcN = len(full_clusters)
    except (KeyError):
        pass

    try:
        fpN = len(set([j for i in full_paths for j in i]))
    except KeyError:
        pass

    logger.info("%s: %s unitigs are created" % (edge,str(fcN+fpN)))
    othercl = len(clusters)-fcN-fpN
    stats.write(edge + "\t" + str(fcN) + "\t" + str(fpN) + "\t" + str(othercl) +"\n")
    stats.close()


def graph_link_unitigs(edge, graph, nx_graph,  bam_cache, link_clusters, link_clusters_src,
                       link_clusters_sink, remove_clusters):
    """
    Second part of the transformation: linkage of all new unitigs created during the first tranforming stage
    """
    logger.debug(f"Linking {edge}")
    link_added = False

    clusters = link_clusters[edge]

    try:
        cl = pd.read_csv("%s/clusters/clusters_%s_%s_%s.csv" % (StRainyArgs().output_intermediate, edge, I, StRainyArgs().AF), keep_default_na = False)
    except(FileNotFoundError):
        pass
    link_unitigs = []

    for phase_clust in set(clusters):
        try:
            if graph.try_get_segment("%s_%s" % (edge, phase_clust)):
                link_unitigs.append(phase_clust)
        except:
            continue

    #for each cluster in the initial unitig
    for cur_clust in link_unitigs:
        #print(f"PROCESSING incoming cluster {cur_clust}")
        cluster_reads = list(cl.loc[cl["Cluster"] == cur_clust, "ReadName"])
        neighbours = {}
        orient = {}

        #get split reads, identify which unitigs they connect and in which orientation
        read_data = bam_cache[edge]
        for read in cluster_reads:
            for next_seg, link_orientation in read_data[read]["Rclip"]:
                try:
                    if len(nx.shortest_path(nx_graph, next_seg, edge)) <= max_hops:
                        neighbours[read] = next_seg
                except (nx.NetworkXNoPath, nx.NodeNotFound):
                    pass

                if link_orientation == "+":
                    orient[next_seg] = ("+", "+")
                else:
                    orient[next_seg] = ("+", "-")

            for next_seg, link_orientation in read_data[read]["Lclip"]:
                try:
                    if len(nx.shortest_path(nx_graph, next_seg, edge)) <= max_hops:
                        neighbours[read] = next_seg
                except (nx.NetworkXNoPath, nx.NodeNotFound):
                    pass

                if link_orientation == "+":
                    orient[next_seg] = ("-", "-")
                else:
                    orient[next_seg] = ("-", "+")

        #print("Neighbors", neighbours)

        #for each "neighbor" (a potential unitig-unitig connection derived from reads)
        for next_seg in set({k for k, v in Counter(neighbours.values()).items() if v >= min_reads_neighbour}):
            #print(f"\tPROCESSING outgoing segment {next_seg}")
            fr_or, to_or = orient[next_seg]
            try:
                cl_n = pd.read_csv("%s/clusters/clusters_%s_%s_%s.csv" % (StRainyArgs().output_intermediate, next_seg,
                                                                          I, StRainyArgs().AF), keep_default_na = False)
            except(FileNotFoundError):
                gfa_ops.add_link(graph, f"{edge}_{cur_clust}", fr_or, next_seg, to_or, 555)
                continue

            #for each neighbor, identify which clusters should be connected
            connecting_reads = []
            for read, read_adj in neighbours.items():
                if read_adj == next_seg:
                    connecting_reads.append(read)
            connected_clusters = cl_n.loc[cl_n["ReadName"].isin(connecting_reads), "Cluster"]
            connected_clusters_thld = list(set([x for x in list(Counter(list(connected_clusters)))
                                            if Counter(list(connected_clusters))[x]  >= min_reads_cluster]))
            #print("Connected clusters", connected_clusters)
            #print("Clusters thld", connected_clusters_thld)

            #make cluster-cluster connections
            link_added = False
            for next_clust in connected_clusters_thld:
                w = Counter(list(connected_clusters))[next_clust]
                try:
                    if graph.try_get_segment(f"{next_seg}_{next_clust}"):
                        gfa_ops.add_link(graph, f"{edge}_{cur_clust}", fr_or, f"{next_seg}_{next_clust}", to_or, w)
                        #print(f"Direct link: {edge}_{cur_clust} to {next_seg}_{next_clust}, {w}")
                        link_added = True
                except(gfapy.NotFoundError):
                    continue

            #in case nothing was connected, either connect it to all starts/ends for the other uinig
            #or, if the partent segment is not deleted, connect to parent segment
            if link_added == False:
                if next_seg in remove_clusters:
                    rewire_clusters = []
                    try:
                        if cur_clust in link_clusters_sink[edge] and cur_clust in link_clusters_src[edge]:
                            if to_or == "+":
                                rewire_clusters = link_clusters_src[next_seg]
                            else:
                                rewire_clusters = link_clusters_sink[next_seg]
                        elif cur_clust in link_clusters_sink[edge]:
                            rewire_clusters = link_clusters_src[next_seg]  #
                        elif cur_clust in link_clusters_src[edge]:
                            rewire_clusters = link_clusters_sink[next_seg]
                    except(KeyError):
                        pass

                    for next_clust in rewire_clusters:
                        try:
                            if graph.try_get_segment(f"{next_seg}_{next_clust}"):
                                #link_added = True
                                gfa_ops.add_link(graph, f"{edge}_{cur_clust}", fr_or, f"{next_seg}_{next_clust}", to_or, 666)
                        except(gfapy.NotFoundError):
                            pass

                else:
                    gfa_ops.add_link(graph, f"{edge}_{cur_clust}", fr_or, next_seg, to_or, 666)
                    #link_added = True
            ###### end block of non-connection


def connect_parental_edges(graph, link_clusters_src, link_clusters_sink, remove_clusters):
    def is_right_tip(seg, sign):
        if graph.segment(seg) is None:
            return False
        if sign == "+":
            return len(graph.segment(seg).dovetails_R) == 0
        else:
            return len(graph.segment(seg).dovetails_L) == 0

    def neg_sign(sign):
        return "+" if sign == "-" else "-"

    #if we keeping the parental segment,
    #connect parent segment with all clusters in the adjecent edge.
    from_index = defaultdict(list)
    to_index = defaultdict(list)
    for link in graph.dovetails:
        from_index[link.from_segment].append(link)
        to_index[link.to_segment].append(link)

    for edge in StRainyArgs().edges:
        if edge in remove_clusters:
            continue

        for link in from_index[edge]:
            to_connect = link_clusters_src[link.to_segment.name] if link.to_orient == "+" else \
                         link_clusters_sink[link.to_segment.name]
            for next_clust in to_connect:
                #logger.debug(str(link).replace(link.to_segment.name, f"{link.to_segment.name}_{next_clust}"))
                candidate = f"{link.to_segment.name}_{next_clust}"
                if is_right_tip(candidate, neg_sign(link.to_orient)):
                    gfa_ops.add_link(graph, link.from_segment.name, link.from_orient,
                             candidate, link.to_orient, 888)

        for link in to_index[edge]:
            to_connect = link_clusters_sink[link.from_segment.name] if link.from_orient == "+" else \
                         link_clusters_src[link.from_segment.name]
            for next_clust in to_connect:
                #logger.debug(str(link).replace(link.from_segment.name, f"{link.from_segment.name}_{next_clust}"))
                candidate = f"{link.from_segment.name}_{next_clust}"
                if is_right_tip(candidate, link.from_orient):
                    gfa_ops.add_link(graph, candidate, link.from_orient,
                             link.to_segment.name, link.to_orient, 888)


def clean_graph(g):
    """
    Remove 0len unitigs, virtual  and self links
    :param g:
    :return:
    """
    for line in g.dovetails:
        if line.from_segment == line.to_segment: #TODO do not self links
            g.rm(line)
        #if g.segment(line.from_segment).virtual == True or g.segment(line.to_segment).virtual == True:
        #    g.rm(line)
    for seq in g.segments:  #TODO do not create o len unitigs
        if len(seq.sequence) == 0:
            seq.sequence = "A"
            #seq.dp = 0
    for path in g.paths:
        g.rm(path)


def transform_main(args):
    init_global_args_storage(args)

    if os.path.isdir(StRainyArgs().log_transform):
        shutil.rmtree(StRainyArgs().log_transform)
    os.mkdir(StRainyArgs().log_transform)
    set_thread_logging(StRainyArgs().log_transform, "transform_root", None)

    stats = open("%s/stats_clusters.txt" % StRainyArgs().output_intermediate, "a")
    stats.write("Edge" + "\t" + "Fill Clusters" + "\t" + "Full Paths Clusters" + "\n")
    stats.close()

    """
    Here we put a hard limit on the number of 16 threads. This is because of an issue in CPython implementation
    of multiprocessing that has a hardcoded contant of max 16 threads that can wait for Lock().
    The issue has been fixed in Python 3.11, and once we transtition to this version, this hard limit
    can be removed.
    https://github.com/katerinakazantseva/stRainy/issues/75
    https://github.com/python/cpython/issues/101225
    """
    HARD_LIMIT = 16
    num_threads = min(StRainyArgs().threads, HARD_LIMIT)
    pool = None
    if StRainyArgs().threads != 1:
        pool = multiprocessing.Pool(num_threads)
    default_manager = multiprocessing.Manager()

    initial_graph = gfapy.Gfa.from_file(StRainyArgs().gfa)
    #Setting up coverage for all unitigs based on bam alignment depth
    logger.info("Re-setting unitigs coverage")
    ref_coverage = {}
    for edge in StRainyArgs().edges:
        edge_cov = pysam.samtools.coverage("-r", edge, StRainyArgs().bam, "--no-header").split()[6]
        initial_graph.try_get_segment(edge).dp = round(float(edge_cov))
        ref_coverage[edge] = round(float(edge_cov))

    logger.info("Loading phased unitigs dictionary")
    try:
        with open(os.path.join(StRainyArgs().output_intermediate, consensus_cache_path), "rb") as f:
            logger.debug(f"searching consensus cache in {os.getcwd()}")
            consensus_dict = pickle.load(f)
    except FileNotFoundError:
        consensus_dict = {}

    flye_consensus = FlyeConsensus(StRainyArgs().bam, StRainyArgs().fa, args.threads, consensus_dict, default_manager)
    consensus_dict = {}

    logger.info("### Create unitigs")
    bam_cache, link_clusters, link_clusters_src, link_clusters_sink, remove_clusters, initial_graph = \
            parallelize_gcu(pool, StRainyArgs().edges, flye_consensus, initial_graph, args)

    # Save phased and reference unitigs' info as a csv
    logger.info('Creating csv file with phased unitigs...')
    write_phased_unitig_csv()
    #logger.info('Done!')
    logger.info('Creating csv file with reference unitigs...')
    store_reference_unitig_info(ref_coverage)
    write_reference_unitig_csv()
    #logger.info('Done!')

    logger.info("### Link unitigs")
    nx_graph = gfa_ops.gfa_to_nx(initial_graph)
    for edge in StRainyArgs().edges:
        graph_link_unitigs(edge, initial_graph, nx_graph, bam_cache, link_clusters, link_clusters_src,
                           link_clusters_sink, remove_clusters)
    connect_parental_edges(initial_graph, link_clusters_src, link_clusters_sink, remove_clusters)

    logger.info("### Remove initial segments")
    for ed in initial_graph.segments:
        if ed.name in remove_clusters:
            logger.debug(f"Removing {ed.name}")
            initial_graph.rm(ed)
    for link in initial_graph.dovetails:
        if link.to_segment in remove_clusters or link.from_segment in remove_clusters:
            initial_graph.rm(link)

    clean_graph(initial_graph)
    out_clusters = os.path.join(StRainyArgs().output_intermediate, "10_fine_clusters.gfa")
    gfapy.Gfa.to_file(initial_graph, out_clusters)

    strainy_utgs = os.path.join(StRainyArgs().output, "strain_unitigs.gfa")
    shutil.copyfile(out_clusters, strainy_utgs)

    phased_graph = gfapy.Gfa.from_file(out_clusters)    #parsing again because gfapy can"t copy
    gfapy.GraphOperations.merge_linear_paths(phased_graph)
    clean_graph(phased_graph)
    out_merged = os.path.join(StRainyArgs().output_intermediate, "20_extended_haplotypes.gfa")
    gfapy.Gfa.to_file(phased_graph, out_merged)

    strainy_final = os.path.join(StRainyArgs().output, "strain_contigs.gfa")
    shutil.copyfile(out_merged, strainy_final)

    logger.info("### Simplify graph")
    smpl.simplify_links(initial_graph)
    gfapy.GraphOperations.merge_linear_paths(initial_graph)
    clean_graph(initial_graph)

    out_simplified = os.path.join(StRainyArgs().output_intermediate, "30_links_simplification.gfa")
    gfapy.Gfa.to_file(initial_graph, out_simplified)
    if args.link_simplify:
        shutil.copyfile(out_merged, strainy_final)

    logger.info("Generating strain report")
    strains_report = os.path.join(StRainyArgs().output, "multiplicity_stats.txt")
    strain_stats_report(StRainyArgs().reference_unitig_info_table_path,
                        StRainyArgs().phased_unitig_info_table_path, open(strains_report, "w"))

    logger.info("Generating strain variant calls")
    strain_utgs_fasta = os.path.join(StRainyArgs().output_intermediate, "strain_utgs.fasta")
    strain_utgs_aln = strain_utgs_fasta + "_ref_aln.bam"
    gfa_to_fasta(out_clusters, strain_utgs_fasta)
    vcf_strain_variants = os.path.join(StRainyArgs().output, "strain_variants.vcf")
    produce_strainy_vcf(StRainyArgs().fa, strain_utgs_fasta, StRainyArgs().threads,
                        strain_utgs_aln, open(vcf_strain_variants, "w"))

    flye_consensus.print_cache_statistics()
    logger.info("### Done!")

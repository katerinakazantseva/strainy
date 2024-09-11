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
from strainy.clustering import build_adj_matrix as matrix
from strainy.clustering import cluster_postprocess as postprocess
from strainy.simplification import simplify_links as smpl
from strainy.graph_operations import gfa_ops
from strainy.graph_operations import asm_graph_ops
from strainy.graph_operations import overlap_graph_ops
from strainy.unitig_statistics import utg_stats
from strainy.flye_consensus import FlyeConsensus
from strainy.clustering import build_data
from strainy.params import *
from strainy.logging import set_thread_logging
from strainy.reports.strainy_stats import strain_stats_report
from strainy.reports.call_variants import produce_strainy_vcf
from strainy.preprocessing import gfa_to_fasta
from strainy.phase import color_bam
from dataclasses import dataclass
logger = logging.getLogger()






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
    """
    Parallelizes the process of unitig creation using worker threads to improve performance,
    especially when handling a large number of graph edges.
    This function manages multithreading to create new unitigs in parallel, and then merges
    the results to update the graph and handle further graph operations.
    """
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
            asm_graph_ops.add_child_edge(op[1], op[2], graph, op[3], op[4], op[5], op[6], flye_consensus, op[7], op[8])
        elif op[0] == 'add_path_edges':
            overlap_graph_ops.add_path_edges(op[1], graph, op[2], op[5], op[6], op[7], op[8], op[9], op[10], op[11], flye_consensus)
        elif op[0] == 'add_path_links':
            asm_graph_ops.add_path_links(graph, op[1], op[2])
        
    return bam_cache, link_clusters, link_clusters_src, link_clusters_sink, set(remove_clusters), graph


def graph_create_unitigs(edge, flye_consensus, bam_cache, link_clusters,
                         link_clusters_src, link_clusters_sink, remove_clusters, graph_ops):
    """
    First stage of the transformation: creation of all new unitigs from clusters obtained during the phasing stage
    Returns:
    None: The function performs operations in-place, modifying GFA graph
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
                    if asm_graph_ops.strong_tail(cluster, cl, ln, data)[0] == True and asm_graph_ops.strong_tail(cluster, cl, ln,data)[1] == True:
                        consensus = flye_consensus.flye_consensus(cluster, edge, cl)
                        # add_child_edge(edge, cluster, graph, cl,consensus["start"], consensus["end"], cons, flye_consensus)
                        graph_ops.append(['add_child_edge', edge, cluster, cl, consensus["start"], consensus["end"], cons, True, True])
                        full_clusters.append(cluster)

                    elif asm_graph_ops.strong_tail(cluster, cl, ln, data)[0] != True:
                        cons[cluster]["Start"] = cons[cluster]["Start"] + start_end_gap+1
                    else:
                        cons[cluster]["End"] = cons[cluster]["End"] - start_end_gap-1
                if clStart < start_end_gap and asm_graph_ops.strong_tail(cluster, cl, ln, data)[0] == True :
                    full_paths_roots.append(cluster)
                if clStop > ln - start_end_gap and asm_graph_ops.strong_tail(cluster, cl, ln, data)[1] == True:
                    full_paths_leafs.append(cluster)

            cluster_distances = postprocess.build_adj_matrix_clusters(edge, cons, cl, flye_consensus, False)
            cluster_distances = matrix.change_w(cluster_distances,0)

            G = overlap_graph_ops.build_paths_graph(cons, full_paths_roots, full_paths_leafs, cluster_distances.copy())

            #full_cl[edge] = full_clusters
            if StRainyArgs().debug:
                overlap_graph_ops.paths_graph_add_vis(edge,
                                    cons,
                                    cl,
                                    full_paths_roots,
                                    full_paths_leafs,
                                    full_clusters,
                                    cluster_distances.copy())

            try:
                full_paths = overlap_graph_ops.find_full_paths(G,full_paths_roots, full_paths_leafs)
            except(ValueError):
                pass

            graph_ops.append(['add_path_edges', edge, cl, data, SNP_pos, ln, full_paths, G,full_paths_roots,
                           full_paths_leafs,full_clusters,cons])
            graph_ops.append(['add_path_links', edge, full_paths])

            othercl = list(set(clusters) - set(full_clusters) - {j for i in full_paths for j in i})
            if len(othercl) > 0:
                G = gfa_ops.from_pandas_adjacency_notinplace(cluster_distances.copy(), create_using = nx.DiGraph)

            close_to_full = []
            othercl_len=[cons[i]['End']-cons[i]['Start'] for i in othercl]
            othercl_sorted=[i[1] for i  in sorted(zip(othercl_len, othercl), reverse=True)]
            removed=[]
            for cluster in othercl_sorted:
                neighbors = nx.all_neighbors(G, cluster)
                A = set(neighbors)
                B = {j for i in full_paths for j in i}
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


            new_cov = asm_graph_ops.change_cov(graph, edge, cons, ln, clusters, othercl, remove_clusters)
            if  new_cov < parental_min_coverage and len(clusters) - len(othercl) != 0 and (len(set(full_clusters))>0 or len(full_paths)>0):
                remove_clusters.add(edge)
            else:
                for cluster in othercl:
                    consensus = flye_consensus.flye_consensus(cluster, edge, cl)
                    # add_child_edge(edge, cluster, graph, cl, cons[cluster]["Start"], cons[cluster]["End"], cons, flye_consensus,insertmain=False)
                    graph_ops.append(['add_child_edge', edge, cluster, cl, cons[cluster]["Start"], cons[cluster]["End"], cons, True, False])
                remove_clusters.add(edge)

            link_clusters[edge] = list(full_clusters) + list(
                set(full_paths_roots).intersection({j for i in full_paths for j in i})) + list(
                set(full_paths_leafs).intersection({j for i in full_paths for j in i}))
            link_clusters_src[edge] = list(full_clusters) + list(
                set(full_paths_roots).intersection({j for i in full_paths for j in i}))
            link_clusters_sink[edge] = list(full_clusters) + list(
                set(full_paths_leafs).intersection({j for i in full_paths for j in i}))

    stats = open("%s/stats_clusters.txt" % StRainyArgs().output_intermediate, "a")
    fcN = 0
    fpN = 0

    try:
        fcN = len(full_clusters)
    except (KeyError):
        pass

    try:
        fpN = len({j for i in full_paths for j in i})
    except KeyError:
        pass

    logger.info("%s: %s unitigs are created" % (edge,str(fcN+fpN)))
    othercl = len(clusters)-fcN-fpN
    stats.write(edge + "\t" + str(fcN) + "\t" + str(fpN) + "\t" + str(othercl) +"\n")
    stats.close()


def graph_link_unitigs(edge, graph, nx_graph,  bam_cache, link_clusters, link_clusters_src,
                       link_clusters_sink, remove_clusters):
    """
    Second part of the transformation: Linking all new unitigs created during the first transformation stage.
    This function links the newly created unitigs (from the first transformation stage) to other unitigs
    based on the read connections between clusters.
    It utilizes BAM data and the graph structure to identify potential connections and adds GFA links between
    unitigs using the `gfa_ops.add_link` function.
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
            connected_clusters_thld = list({x for x in list(Counter(list(connected_clusters)))
                                            if Counter(list(connected_clusters))[x]  >= min_reads_cluster})

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
    """
    Connects the parental segments to the newly created clusters in adjacent edges.
    This function ensures that the original (parental) segments are connected to the newly created clusters
    during the second transformation stage. It identifies segments that should remain in the graph
    (i.e., not in `remove_clusters`) and links them to their neighboring clusters using GFA links.
    The function handles the connections based on segment orientation and
    the source/sink clusters associated with each edge.
    """
    def is_right_tip(seg, sign):
        if graph.segment(seg) is None:
            return False
        if sign == "+":
            return len(graph.segment(seg).dovetails_R) == 0

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
                candidate = f"{link.to_segment.name}_{next_clust}"
                if is_right_tip(candidate, neg_sign(link.to_orient)):
                    gfa_ops.add_link(graph, link.from_segment.name, link.from_orient,
                             candidate, link.to_orient, 888)

        for link in to_index[edge]:
            to_connect = link_clusters_sink[link.from_segment.name] if link.from_orient == "+" else \
                         link_clusters_src[link.from_segment.name]
            for next_clust in to_connect:
                candidate = f"{link.from_segment.name}_{next_clust}"
                if is_right_tip(candidate, link.from_orient):
                    gfa_ops.add_link(graph, candidate, link.from_orient,
                             link.to_segment.name, link.to_orient, 888)



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
    utg_stats.write_phased_unitig_csv()
    logger.info('Creating csv file with reference unitigs...')
    utg_stats.store_reference_unitig_info(ref_coverage)
    utg_stats.write_reference_unitig_csv()


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

    initial_graph=gfa_ops.clean_graph(initial_graph)
    out_clusters = os.path.join(StRainyArgs().output_intermediate, "10_fine_clusters.gfa")
    gfapy.Gfa.to_file(initial_graph, out_clusters)

    strainy_utgs = os.path.join(StRainyArgs().output, "strain_unitigs.gfa")
    shutil.copyfile(out_clusters, strainy_utgs)

    phased_graph = gfapy.Gfa.from_file(out_clusters)    #parsing again because gfapy can"t copy

    segs_unmerged=phased_graph.segment_names
    gfapy.GraphOperations.merge_linear_paths(phased_graph)
    phased_graph=gfa_ops.clean_graph(phased_graph)
    segs_merged = phased_graph.segment_names

    out_merged = os.path.join(StRainyArgs().output_intermediate, "20_extended_haplotypes.gfa")
    gfapy.Gfa.to_file(phased_graph, out_merged)

    strainy_final = os.path.join(StRainyArgs().output, "strain_contigs.gfa")
    shutil.copyfile(out_merged, strainy_final)

    logger.info("### Simplify graph")
    smpl.simplify_links(initial_graph)
    gfapy.GraphOperations.merge_linear_paths(initial_graph)
    initial_graph=gfa_ops.clean_graph(initial_graph)

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

    logger.info("Update clusters (csv) and colored BAM")
    merged_clusters={}
    AF = StRainyArgs().AF
    for seg in [i for i in segs_unmerged if i not in segs_merged]:
        seg_merged = [k for k in segs_merged if re.search(seg, k) != None][0]
        merged_clusters[seg] = seg_merged

    for edge in StRainyArgs().edges:
        try:
            cl = pd.read_csv("%s/clusters/clusters_%s_%s_%s.csv" % (StRainyArgs().output_intermediate, edge, I, AF),
                         keep_default_na=False)
            clusters = sorted(set(cl['Cluster']))
            for cluster in clusters:
                seg=str(edge)+"_"+str(cluster)
                if  seg in merged_clusters.keys():
                    cl.loc[cl['Cluster'] == cluster, 'Cluster'] = merged_clusters[seg]
            cl.to_csv("%s/clusters/clusters_%s_%s_%s_MERGED.csv" % (StRainyArgs().output_intermediate, edge, I, AF))
        except(FileNotFoundError): pass
        os.makedirs("%s/bam/merged/" % StRainyArgs().output_intermediate, exist_ok=True)
    color_bam(StRainyArgs().edges, transfrom_stage=True)
    flye_consensus.print_cache_statistics()
    logger.info("### Done!")
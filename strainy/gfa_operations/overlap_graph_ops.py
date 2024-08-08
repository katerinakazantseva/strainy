import networkx as nx
import pygraphviz as gv
from collections import Counter, deque, defaultdict
import logging
from strainy.gfa_operations import gfa_ops
from strainy.gfa_operations import asm_graph_ops
from strainy.params import *


logger = logging.getLogger()

"""
This contains functions for operation with overlap graph:
1. build_paths_graph: creates overlap graph  #TODO rename to overlap
2. find_full_paths: finds full paths in overlap graph
3.remove_nested(G, cons): removes nested clusters (
4.add_path_edges : calc cluster boundaries and creates unitigs using asm.add_child_edge

"""


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
    #todo move it to parental function
    G, full_paths_roots, full_paths_leafs = \
        remove_leaf_root_subnodes(G,full_paths_roots,full_paths_leafs)
    G = remove_nested(G, cons)
    G = remove_transitive(G)
    return G


def remove_transitive(G):
    path_remove = []
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
    return G


def remove_leaf_root_subnodes(G,full_paths_roots,full_paths_leafs):
    node_remove = []
    for node in full_paths_leafs+full_paths_roots:
        if node in full_paths_leafs:
            neighbors = list(full_paths_leafs)
        else:
            neighbors = list(full_paths_roots)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G, node, neighbor, cutoff = 2):
                if len(n_path) == 2:
                    node_remove.append(neighbor)
    for node in node_remove:
        try:
            G.remove_node(node)
            logger.debug("REMOVE " + str(node))
            full_paths_roots.remove(node)
            full_paths_leafs.remove(node)
        except:
            continue
    return (G,full_paths_roots,full_paths_leafs)


def remove_bubbles(graph, source_nodes):
    for node in source_nodes:
        neighbors = list(source_nodes)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(graph, node, neighbor, cutoff = 2):
                if len(n_path) == 2:
                    node_remove.append(neighbor)




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

    return paths



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
    return G

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
            asm_graph_ops.add_child_edge(edge, path_cluster, g,  cl, cut_l[path_cluster], cut_r[path_cluster], cons, flye_consensus)
        else:
            for i in range(0,len(full_paths)):
                if path_cluster in full_paths[i]:
                    upd_path = full_paths[i]
                    upd_path.remove(path_cluster)
                    full_paths[i] = upd_path
            G.remove_node(path_cluster)

    return(path_cl)


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



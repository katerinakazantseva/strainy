import csv
import networkx as nx
from build_adj_matrix import *
from cluster_postprocess import *
import pygraphviz as gv
import re
import gfapy
from collections import Counter, deque
from build_data  import *
from params import *
import numpy as np
from simplify_links import *

g = gfapy.Gfa.from_file(gfa)
stats = open('output/stats_clusters.txt', 'a')
stats.write("Edge" + "\t" + "Fill Clusters" + "\t" + "Full Paths Clusters" + "\n")
stats.close()



def add_child_edge(edge, clN, g, cl, SNP_pos, data, left, righ,cons):
    seq=[]
    cl_consensuns = cluster_consensuns(cl, clN, SNP_pos, data, cons)
    try:
        i=g.try_get_segment(edge)
        seq = i.sequence
        seq = list(seq)

        for key, val in cl_consensuns[clN].items():
            try:
                seq[int(key)-1] = val
            except (ValueError, KeyError):
                continue
        seq = ''.join(seq)
        seq = seq[left:righ+1]
    except(gfapy.NotFoundError):
        pass

    if len(seq)==0:
        remove_zeroes.append("S\t%s_%s\t*" % (edge, clN))
    if len(seq)>0:
        g.add_line("S\t%s_%s\t*" % (edge, clN))
        i = g.try_get_segment("%s_%s" % (edge, clN))
        new_line = i
        new_line.name = str(edge) + "_" + str(clN)
        new_line.sid = str(edge) + "_" + str(clN)
        new_line.sequence = seq
        new_line.dp = cons[clN]["Cov"]  # coverage
        print("edge added:" + str(new_line.name))

'''

def remove_link(fr,fr_or, to, to_or):
    res=False
    for i in g.dovetails:
        if (i.from_segment == fr and i.to_segment == to):
            g.rm(i)
            res=True
            print("remove line: " + str(i.from_segment) + str(i.from_orient) + " to " + str(i.to_segment) + str(i.to_orient))
    return (res)'''



def build_paths_graph(SNP_pos, cl, cons,full_clusters, data,ln, full_paths_roots, full_paths_leafs):

    M = build_adj_matrix_clusters(cons, SNP_pos, cl, False)
    M = change_w(M, 1)
    G = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    try:
        G.remove_node(0)
    except:
        pass

    """
    for e in G.edges():
        first_cl=e[0]
        second_cl=e[1]
        intersect = set(range(cons[first_cl]["Start"], cons[first_cl]["Stop"])).intersection(
            set(range(cons[second_cl]["Start"], cons[second_cl]["Stop"])))
        G[e[0]][e[1]]['weight'] = len(intersect)"""
    path_remove = []
    node_remove=[]

    for node in full_paths_leafs:
        neighbors = list(full_paths_leafs)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G, node, neighbor):
                if len(n_path) == 2:
                    node_remove.append(neighbor)

    for node in full_paths_roots:
        neighbors = list(full_paths_roots)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G,  neighbor,node):
                if len(n_path) == 2:
                    node_remove.append(neighbor)


    '''
    for full_cluster in full_clusters:
        if strong_tail(full_cluster, cl, ln, "root", data)==True and strong_tail(full_cluster, cl, ln, "leaf", data)==True:
            G.remove_node(full_cluster)
            full_paths_roots.remove(full_cluster)
            full_paths_leafs.remove(full_cluster)
            #print("------------------")'''

    G=remove_nested(G, cons)
    #G = remove_not_strong_tails(G, cl, cons, data, ln, full_paths_roots, full_paths_leafs)


    for node in node_remove:
        try:
            G.remove_node(node)
            print("REMOVE "+str(node))
            full_paths_roots.remove(node)
            full_paths_leafs.remove(node)
        except:
            continue


    for node in G.nodes():
        neighbors = nx.all_neighbors(G, node)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G, node, neighbor):
                if len(n_path) == 3:
                    path_remove.append(n_path)


    for n_path in path_remove:
        #print(n_path)
        try:
            #G.remove_edge(n_path[0], n_path[2])
            G.remove_edge(n_path[0], n_path[1])
            #G.remove_edge(n_path[1], n_path[2])


        except:
            continue

    return (G)

def remove_nested(G, cons):
    nodes=list(G.nodes())
    #print(G.nodes)
    for node in nodes:
        try:
            neighbors = nx.all_neighbors(G, node)
            for neighbor in list(neighbors):
                if cons[node]["Start"]<cons[neighbor]["Start"] and cons[node]["Stop"]>cons[neighbor]["Stop"]:
                    try:
                        G.remove_node(neighbor)
                        print("REMOVE NESTED"+str(neighbor))
                    except:
                        continue
        except:
            continue
    #print(G.nodes)
    return (G)


def paths_graph_add_vis(edge,cons, SNP_pos, cl,full_paths_roots,full_paths_leafs,full_clusters):
    M = build_adj_matrix_clusters(cons, SNP_pos, cl, False)
    M = change_w(M, 1)
    G_vis = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    cl_removed = []
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    try:
        G_vis.remove_node(0)
    except:
        pass
    path_remove = []

    G_vis=remove_nested(G_vis, cons)
    for node in G_vis.nodes():
        neighbors = nx.all_neighbors(G_vis, node)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G_vis, node, neighbor):
                if len(n_path) == 3:
                    path_remove.append(n_path)



    for e in G_vis.edges():
        first_cl = e[0]
        second_cl = e[1]
        intersect = set(range(cons[first_cl]["Start"], cons[first_cl]["Stop"])).intersection(
            set(range(cons[second_cl]["Start"], cons[second_cl]["Stop"])))
        G_vis[e[0]][e[1]]['weight'] = len(intersect)


    G_vis.add_node("Src",style='filled',fillcolor='gray',shape='square')
    G_vis.add_node("Sink",style='filled',fillcolor='gray',shape='square')
    for i in full_paths_roots:
        G_vis.add_edge("Src", i)
    for i in full_paths_leafs:
        G_vis.add_edge(i, "Sink")
    for i in full_clusters:
        G_vis.add_edge("Src", i)
        G_vis.add_edge(i, "Sink")

    test = nx.nx_agraph.to_agraph(G_vis)

    test = str(test)
    #test = test.replace('weight', 'label')
    test = gv.AGraph(test)
    test.layout(prog="neato")
    test.draw("output/graphs/full_paths_cluster_GV_graph_%s.png" % (edge))
    G_vis.remove_node("Src")
    G_vis.remove_node("Sink")
    return(cl_removed)

def find_full_paths(G, paths_roots, paths_leafs):
    paths = []
    for root in paths_roots:
        paths_nx = nx.algorithms.all_simple_paths(G, root, paths_leafs)
        for path in list(paths_nx):
            paths.append(path)
    return (paths)


def add_link(fr, fr_or, to, to_or,w):
    link = 'L	%s	%s	%s	%s	0M	ex:i:%s' % (fr, fr_or, to, to_or, w)
    #print(link)
    try:
        g.add_line(link)
        print("link added from %s %s to %s %s" % (fr, fr_or, to, to_or))
    except(gfapy.NotUniqueError): pass
        #print("dd")


def add_path_links(edge, paths, G):
    #str = 'L	first_edge	+	second_edge	+	0M	RC:s:in'
    for path in paths:
        for i in range(0, len(path) - 1):
                try:
                    #w=G[path[i]][path[i+1]]['weight']
                    str='L	first_edge	+	second_edge	+	0M	ix:i:%s' % 1
                    #print("path added between clusters" + path[i] + " and " + path[i + 1])
                    g.add_line(str.replace('first_edge', "%s_%s" % (edge, path[i])).replace('second_edge',
                                                                                            "%s_%s" % (
                                                                               edge, path[i + 1])))
                except(gfapy.error.NotUniqueError, KeyError):
                    continue





def add_path_edges ( edge,g,cl, data, SNP_pos, ln, paths, G,paths_roots,paths_leafs, cons):

    path_cl = []

    for path in paths[edge]:
        for member in path:
            path_cl.append(member)
    cut_l_unsorted = {}
    cut_r_unsorted = {}




    for path_cluster in set(path_cl):
        cut_l_unsorted[path_cluster] = None
        cut_r_unsorted[path_cluster] = None
        if path_cluster in paths_roots and cons[path_cluster]["Start"]<3 :

            cut_l_unsorted[path_cluster] = 0
            #cut_r[path_cluster] = None
        if path_cluster in paths_leafs:

            cut_r_unsorted[path_cluster] = ln - 1
            #cut_l[path_cluster] = None
        #else:
            #cut_l[path_cluster] = None
            #cut_r[path_cluster] = None
    stop_pos={}
    for i in cut_r_unsorted.keys():
        stop_pos[i]=cons[i]["Stop"]


    order_by_stop_pos = list(dict(sorted(stop_pos.items(), key=lambda item: item[1])).keys())

    cut_l = {}
    cut_r = {}
    for i in order_by_stop_pos:
        cut_l[i] = cut_l_unsorted[i]
        cut_r[i] = cut_r_unsorted[i]


    while None in cut_l.values():
        for member in cut_l.keys():
            if cut_l[member] != None and (cut_r[member] == None or member in paths_leafs):
                #print(member)
                Q = deque()
                L = []
                R = []
                for path in paths[edge]:

                    try:
                        L.append(path[path.index(member) + 1])

                        Q.append(path[path.index(member) + 1])
                    except (ValueError, IndexError):
                        continue
                visited = []
                while Q:

                    n = Q.pop()

                    visited.append(n)
                    if n in L:
                        for path in paths[edge]:
                            try:
                                if path.index(n)>0:
                                    #R.append(path[path.index(n) - 1])
                                    if path[path.index(n) - 1] not in visited:
                                        #print("add R" + path[path.index(n) - 1])
                                        R.append(path[path.index(n) - 1])
                                        Q.append(path[path.index(n) - 1])
                            except (ValueError, IndexError):
                                continue
                    else:
                        for path in paths[edge]:
                            try:
                                #L.append(path[path.index(n) + 1])
                                if path[path.index(n) + 1] not in visited:
                                    L.append(path[path.index(n) + 1])
                                    #print("add L" + path[path.index(n) + 1])
                                    Q.append(path[path.index(n) + 1])

                            except (ValueError, IndexError):
                                continue
                l_borders = []
                r_borders = []
                for i in L:
                    l_borders.append(int(cons[i]["Start"]))
                for i in R:
                    r_borders.append(int(cons[i]["Stop"]))
                if member in paths_leafs:
                    border=cut_r[member]
                else: border = max(l_borders) + (min(r_borders) - max(l_borders)) // 2

                L=list(set(L))
                R = list(set(R))
                #print(L)
                #print(R)

                for i in L:
                    cut_l[i] = border

                for i in R:
                    cut_r[i] = border



    for path_cluster in set(path_cl):
        if cut_l[path_cluster]!=cut_r[path_cluster]:
            add_child_edge(edge, path_cluster, g,  cl, SNP_pos, data, cut_l[path_cluster], cut_r[path_cluster], cons)
        else:
            for i in range(0,len(paths[edge])):
                if path_cluster in paths[edge][i]:
                    upd_path=paths[edge][i]
                    upd_path.remove(path_cluster)
                    paths[edge][i]=upd_path
            G.remove_node(path_cluster)
    return(path_cl)



def change_cov(g,edge,cons,ln,clusters,othercl):
    cov=0
    len_cl=[]
    for i in othercl:
        cov=cov+cons[i]["Cov"]*(cons[i]["Stop"]-cons[i]["Start"])
        for i in range(cons[i]["Start"],cons[i]["Stop"]):
            len_cl.append(i)
    if (len(set(len_cl))/ln)<0.7 and len(clusters)-len(othercl)!=0:
        remove_clusters.append(edge)
    cov=cov/ln
    i = g.try_get_segment(edge)
    i.dp =round(cov)
    return(cov)

def change_sec(g,edge, othercl, cl,SNP_pos, data, cut=True):
    temp={}
    other_cl=cl
    for cluster in othercl:
        other_cl.loc[cl['Cluster']==cluster, "Cluster"] = "OTHER_%s" %edge
    cl_consensuns = cluster_consensuns(other_cl, "OTHER_%s" %edge, SNP_pos, data, temp)

    i = g.try_get_segment(edge)
    seq = i.sequence
    seq = list(seq)
    for key, val in cl_consensuns["OTHER_%s" %edge].items():
        try:
            #seq[int(key)] = val
            seq[int(key) - 1] = val
            #print("changed")
        except (ValueError):
            continue


    seq = ''.join(seq)
    if cut==True:
        #print("CUTTING")
        #print(cl_consensuns["OTHER_%s" % edge]["Start"])
        seq = seq[cl_consensuns["OTHER_%s" % edge]["Start"]:cl_consensuns["OTHER_%s" % edge]["Stop"] + 1]
    i.sequence=seq



'''
def to_neighbours(g,edge,orient):
    to_ng=[]
    for i in g.dovetails:
            if i.from_segment.name==edge and i.from_orient=='+':
                nei=[i.to_segment.name,i.to_orient]
                to_ng.append(nei)
                #print(i)
            if i.to_segment.name==edge and i.to_orient=='-':
                nei=[i.from_segment.name,i.from_orient]
                to_ng.append(nei)
                #print(i)
    return (to_ng)


def from_neighbours(g,edge, orient):
    from_ng=[]
    for i in g.dovetails:

            if i.to_segment.name==edge and i.to_orient=='+':
                nei = [i.from_segment.name, i.from_orient]
                from_ng.append(nei)

            if i.from_segment.name == edge and (i.from_orient=='-'):
                nei = [i.to_segment.name, i.to_orient]
                from_ng.append(nei)
    return (from_ng) '''


def cut(edge):
    fr=1
    to=1
    for i in to_neighbours(g,edge,'+'):
        if len(from_neighbours(g,i[0],i[1]))==1:
            fr=0
    for i in to_neighbours(g,edge,'-'):
        if len(from_neighbours(g,i[0],i[1]))==1:
            fr=0
    for i in from_neighbours(g, edge,'+'):
        if len(to_neighbours(g, i[0],i[1])) == 1:
            to=0
    for i in from_neighbours(g, edge,'-'):
        if len(to_neighbours(g, i[0],i[1])) == 1:
            to=0
    if fr ==0 or to==0:
        res=True
    else: res=False
    return res



def strong_tail(cluster, cl, ln, data):
    count_start = 10000
    count_stop = 10000
    res=[False,False]
    reads = list(cl.loc[cl['Cluster'] == cluster, 'ReadName'])
    for read in reads:
        if data[read]["Start"] < 5:
            if count_start == 10000:
                count_start = 0
            count_start=count_start+1
        if data[read]["Stop"] == ln - 5:
            if count_stop == 10000:
                count_stop = 0
            count_stop = count_stop + 1
    if count_start>2:
        res[0] = True
    if count_stop>2:
        res[1] = True
    return (res)

'''
def remove_not_strong_tails(G, cl,cons, data,ln, full_paths_roots, full_paths_leafs):
    path_remove=[]
    for cluster in full_paths_roots:
        if strong_tail(cluster, cl, ln, "root", data) == True:
            neighbors = nx.all_neighbors(G, cluster)
            for neighbor in list(neighbors):
                for n_path in nx.algorithms.all_simple_paths(G, neighbor, cluster):
                    if len(n_path) == 2:
                        path_remove.append(n_path)

        else:
            #print(cluster)
            full_paths_roots.remove(cluster)
            #print(cons[cluster]["Start"])
            if strong_tail(cluster, cl, ln, "root", data) == 'UnsufStart':
                cons[cluster]["Start"] = cons[cluster]["Start"] + 2
            if strong_tail(cluster, cl, ln, "root", data) == 'UnsufStop':
                cons[cluster]["Stop"] = cons[cluster]["Stop"]-2
            #cons[cluster]["Start"]=2
    for cluster in full_paths_leafs:
        if strong_tail(cluster, cl, ln, "leaf", data) == True:
            neighbors = nx.all_neighbors(G, cluster)
            for neighbor in list(neighbors):
                for n_path in nx.algorithms.all_simple_paths(G, cluster, neighbor):
                    if len(n_path) == 2:
                        path_remove.append(n_path)
        else:

            full_paths_leafs.remove(cluster)
            #print(cons[cluster]["Stop"])
            #с какой стороны бедный??
            if strong_tail(cluster, cl, ln, "leaf", data) == 'UnsufStart':
                cons[cluster]["Start"] = cons[cluster]["Start"] + 2
            if strong_tail(cluster, cl, ln, "leaf", data) == 'UnsufStop':
                cons[cluster]["Stop"] = cons[cluster]["Stop"]-2
            #print(cons[cluster]["Stop"])
    for n_path in path_remove:
        #print(n_path)
        G.remove_edge(n_path[0], n_path[1])


    return (G)
'''

def gfa_to_nx(g):
    G = nx.Graph()
    for i in g.segment_names:
        G.add_node(i)
    for i in g.dovetails:
        G.add_edge(i.from_segment.name, i.to_segment.name)
    return(G)

full_cl = {}
full_paths = {}
full_paths_leafs_roots = {}
paths = {}
full_path_clusters = {}
subunits_borderline = {}
connected_subunits = {}
link_clusters = {}
link_clusters_src = {}
link_clusters_sink = {}
remove_clusters = []
remove_zeroes = []
all_data={}

def graph_create_unitigs(i):
    edge = edges[i]
    print(edge)
    full_paths_roots = []
    full_paths_leafs = []
    full_clusters = []
    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
        SNP_pos = read_snp(snp, edge,bam, AF)
        # Save
        try:
            data=all_data[edge]
        except(KeyError, FileNotFoundError):
            data = read_bam(bam, edge, SNP_pos, clipp, min_mapping_quality, min_al_len, de_max)
            all_data[edge]=data

        ln = int(pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4])
        clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))

        try:
            clusters.remove(0)
        except:
            pass
        cons = build_data_cons(cl, SNP_pos, data)
        #print(clusters)
        if len(clusters) == 1:
            for cluster in clusters:
                clStart=cons[cluster]["Start"]
                clStop = cons[cluster]["Stop"]
                if clStart < 3 and clStop == ln - 1:

                    full_paths_roots.append(cluster)
                    full_paths_leafs.append(cluster)

                add_child_edge(edge, cluster, g, cl, SNP_pos, data, 0, ln - 1, cons)
            link_clusters[edge] = list(clusters)
            link_clusters_sink[edge] = list(clusters)
            link_clusters_src[edge] = list(clusters)
            remove_clusters.append(edge)

        if len(clusters)>1:
            for cluster in clusters:

                clStart=cons[cluster]["Start"]
                clStop = cons[cluster]["Stop"]
                if clStart < 5 and clStop > ln - 5:
                    #full_paths_roots.append(cluster)
                    #full_paths_leafs.append(cluster)
                    if strong_tail(cluster, cl, ln, data)[0] == True and strong_tail(cluster, cl, ln,
                                                                                         data)[1] == True:
                        add_child_edge(edge, cluster, g, cl, SNP_pos, data, 0, ln - 1, cons)
                        full_clusters.append(cluster)

                    elif strong_tail(cluster, cl, ln, data)[0] != True:
                        cons[cluster]["Start"] = cons[cluster]["Start"] + 10
                    else:
                        cons[cluster]["Stop"] = cons[cluster]["Stop"] - 10
                if clStart < 5 and strong_tail(cluster, cl, ln, data)[0] == True :
                    full_paths_roots.append(cluster)

                if clStop > ln - 5 and strong_tail(cluster, cl, ln, data)[1] == True:
                    full_paths_leafs.append(cluster)

            #G = build_paths_graph(SNP_pos, cl, cons, full_clusters)
            G=build_paths_graph(SNP_pos, cl, cons, full_clusters, data, ln, full_paths_roots, full_paths_leafs)
            #G = remove_not_strong_tails(G, cl, cons,data,ln,full_paths_roots,full_paths_leafs)


            full_cl[edge] = full_clusters
            cl_removed=paths_graph_add_vis(edge,cons, SNP_pos,cl,full_paths_roots, full_paths_leafs,full_clusters)
            try:
                full_paths[edge] = find_full_paths(G,full_paths_roots, full_paths_leafs)
                #print(full_paths)
            except(ValueError):
                pass

            add_path_edges(edge,g,cl, data, SNP_pos, ln,full_paths, G,full_paths_roots, full_paths_leafs,cons)
            #print("paths clusters added")
            #print(G)
            add_path_links(edge, full_paths[edge], G)
            #print("paths links added")

            #print("CREATED clusters")



            othercl=list(set(clusters)-set(full_clusters)-set([j for i in full_paths[edge] for j in i])-set(cl_removed))

            close_to_full=[]
            #print(othercl)
            for cluster in othercl.copy():
                #print(cluster)
                M = build_adj_matrix_clusters(cons, SNP_pos, cl, False)
                M = change_w(M, 1)
                G = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
                neighbors = nx.all_neighbors(G, cluster)
                A=set(neighbors)
                B=set([j for i in full_paths[edge] for j in i])
                #print(A)
                #print(B)
                #print(A.intersection(B))
                if len(set(list(neighbors)).intersection(set(full_clusters)))>0 or len(A.intersection(B))>0:
                    othercl.remove(cluster)
                    close_to_full.append(cluster)
                    print("REMOVE "+str(cluster))

            #if len(othercl)>0:
                #othercl=othercl+close_to_full


            #print(othercl)

            new_cov=change_cov(g,edge,cons,ln,clusters,othercl)
            if new_cov<6 and len(clusters)-len(othercl)!=0: #PARAMETER
            #if len(clusters) - len(othercl) != 0:
                remove_clusters.append(edge)
            else:
                #othercl=othercl+close_to_full
                change_sec(g, edge, othercl, cl, SNP_pos, data,True)

            link_clusters[edge] = list(full_clusters) + list(
                set(full_paths_roots).intersection(set([j for i in full_paths[edge] for j in i]))) + list(
                set(full_paths_leafs).intersection(set([j for i in full_paths[edge] for j in i])))
            link_clusters_src[edge] = list(full_clusters) + list(
                set(full_paths_roots).intersection(set([j for i in full_paths[edge] for j in i])))
            link_clusters_sink[edge] = list(full_clusters) + list(
                set(full_paths_leafs).intersection(set([j for i in full_paths[edge] for j in i])))
        else:
            change_sec(g, edge, [clusters[0]], cl, SNP_pos, data, False)
    except(FileNotFoundError, IndexError):
        print("NO CLUSTERS")
        cov = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[6]
        i = g.try_get_segment(edge)
        print(cov)
        i.dp = round(float(cov))
        pass
        clusters = []
    stats = open('output/stats_clusters.txt', 'a')
    fcN = 0
    fpN = 0

    try:
        fcN = len(full_cl[edge])
    except(KeyError):
        pass
    try:

        fpN = len(set([j for i in full_paths[edge] for j in i]))
    except(KeyError,UnboundLocalError):
        pass

    othercl=len(clusters)-fcN-fpN
    stats.write(edge + "\t" + str(fcN) + "\t" + str(fpN) + "\t" + str(othercl) +"\n")
    stats.close()



def graph_link_unitigs(i,G):
    print("CREATING NEW LINKS")
    edge = edges[i]
    print(edge)
    link_added = False

    clusters=[]
    try:
        clusters = link_clusters[edge]
    except(KeyError):
        pass
    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
    except(FileNotFoundError):
        pass
    link_unitigs=[]


    for clN in set(clusters):
        try:
            if g.try_get_segment("%s_%s" % (edge, clN)):
                link_unitigs.append(clN)
        except:
            continue

    print(link_unitigs)

    #for clN in set(clusters):
    for clN in link_unitigs:
        print()
        print("START")
        print("%s_%s" % (edge, clN))
        reads = list(cl.loc[cl['Cluster'] == clN, 'ReadName'])
        neighbours={}
        orient={}



        if 1==1:
            data=all_data[edge]
            for read in reads:


                    for n, v in data[read]["Rclip"].items():
                        #print(data[read]["Rclip"])
                        try:
                            if len(nx.shortest_path(G,n,edge))<=10:
                                #neighbours[read]=n
                                #orient[n]=['+','+']
                                neighbours[read]=n
                                #neighbours.append(n)
                                #orient[n]=[v[1],v[0]]
                                if v[0]=='+' and v[1]=='+':
                                    orient[n] = ['+', '+']
                                elif v[0] == '-' and v[1] == '-':
                                    orient[n] = ['+', '+']
                                else:
                                    #orient[n]=[v[0],v[1]]
                                    orient[n] = ['+', '-']
                        except(nx.NetworkXNoPath):
                            if v[0] == '+' and v[1] == '+':
                                orient[n] = ['+', '+']
                            elif v[0] == '-' and v[1] == '-':
                                orient[n] = ['+', '+']
                            else:
                                # orient[n]=[v[0],v[1]]
                                orient[n] = ['+', '-']


                    ##for n in data[read]["Lclip"]:
                    for n, v in data[read]["Lclip"].items():
                        #print(data[read]["Lclip"])
                        try:
                            if len(nx.shortest_path(G, n, edge)) <= 10:
                                #neighbours[read]=n
                                #orient[n]=['-','-']
                                neighbours[read]=n
                                #neighbours.append(n)
                                #orient[n] = [v[1], v[0]]

                                if v[0]=='+' and v[1]=='+':
                                    orient[n] = ['-', '-']
                                elif v[0] == '-' and v[1] == '-':
                                    orient[n] = ['-', '-']
                                else:
                                    #orient[n]=[v[1], v[0]]
                                    orient[n] = ['-', '+']
                        except(nx.NetworkXNoPath):
                            if v[0] == '+' and v[1] == '+':
                                orient[n] = ['-', '-']
                            elif v[0] == '-' and v[1] == '-':
                                orient[n] = ['-', '-']
                            else:
                                # orient[n]=[v[1], v[0]]
                                orient[n] = ['-', '+']
                            #pass


        #print(neighbours)
        #print(set(neighbours))
        #print(set(neighbours.values()))

        for n in set({k for k, v in Counter(neighbours.values()).items() if v > 0}):
        #for n in set({k for k, v in Counter(neighbours).items() if v > 4}):
            #print(n)
            #print(orient[n])
            fr_or=orient[n][0]
            to_or=orient[n][1]
            w=1
            try:
                cl_n = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (n, I, AF), keep_default_na=False)
            except(FileNotFoundError):
                add_link("%s_%s" % (edge, clN), fr_or,n, to_or,w)
                continue
            reads = []
            for k, v in neighbours.items():
                if v == n:
                    reads.append(k)
            #print(reads)
            n_cl = cl_n.loc[cl_n['ReadName'].isin(reads), 'Cluster']

            #n_cl_set = list(
                #set([x for x in list(n_cl) if Counter(list(n_cl))[x] / sum(Counter(list(n_cl)).values()) >= 0.001]))
            #print(Counter(list(n_cl)))

            n_cl_set = list(set([x for x in list(Counter(list(n_cl))) if Counter(list(n_cl))[x]  >= 2]))
            print(Counter(list(n_cl)))
            print("nei "+str(n) +"cl "+str(n_cl_set))
            '''
            if len(n_cl_set)==0:
                try:

                    #add_link("%s_%s" % (edge, clN), fr_or, n, to_or)
                    #когда мы не находим соседей, соединяем со всеми
                    if n in remove_clusters:
                        try:
                            if clN in link_clusters_sink[edge]:
                                n_cl_set = link_clusters_src[n] #
                            if clN in link_clusters_src[edge]:
                                n_cl_set = link_clusters_sink[n]
                        except(KeyError):
                            pass
                        #for n_cluster in link_clusters_src[n]:

                    else:
                        add_link("%s_%s" % (edge, clN), fr_or, n, to_or,w)
                        #n_cl_set=link_clusters_src[n]
                except (KeyError):
                    continue
            '''

            link_added=False

            for i in n_cl_set:
                w=Counter(list(n_cl))[i]
                try:
                    if g.try_get_segment("%s_%s" % (n, i)):
                        link_added=True
                        add_link("%s_%s" % (edge, clN), fr_or, "%s_%s" % (n, i), to_or,w)
                except(gfapy.NotFoundError):
                    continue

            if link_added==False:
                if n in remove_clusters:
                    #n_cl_set = link_clusters_src[n]
                    try:
                        if clN in link_clusters_sink[edge]:
                            n_cl_set = link_clusters_src[n]  #
                        if clN in link_clusters_src[edge]:
                            n_cl_set = link_clusters_sink[n]
                    except(KeyError):
                        pass
                    # for n_cluster in link_clusters_src[n]:

                else:
                    add_link("%s_%s" % (edge, clN), fr_or, n, to_or,w)
                    link_added = True
                    #n_cl_set = link_clusters_src[n]
            for i in n_cl_set:
                try:
                    if g.try_get_segment("%s_%s" % (n, i)):
                        link_added=True
                        add_link("%s_%s" % (edge, clN), fr_or, "%s_%s" % (n, i), to_or,w)
                except(gfapy.NotFoundError):
                    continue




    if link_added==False or edge not in remove_clusters:
        print("restore links")
        for d in g.dovetails:
            repl=[]
            if d.from_segment==edge:
                if d.to_orient=='+':
                    try:
                        for i in link_clusters_src[d.to_segment.name]:
                            repl.append(i)
                    except(KeyError):
                        pass
                if d.to_orient == '-':
                    try:
                        for i in link_clusters_sink[d.to_segment.name]:
                            repl.append(i)
                    except(KeyError):
                        pass
                for i in repl:
                    print(str(d).replace(d.to_segment.name,'%s_%s' % (d.to_segment.name,i)))
                    try:
                        g.add_line(str(d).replace(d.to_segment.name,'%s_%s' % (d.to_segment.name,i)))
                    except(gfapy.error.NotUniqueError):
                        pass
            if d.to_segment==edge:
                print(d.from_segment.name)
                if d.from_orient == '+':
                    try:
                        for i in link_clusters_sink[d.from_segment.name]:
                            repl.append(i)
                    except(KeyError):
                        pass
                if d.from_orient == '-':
                    try:
                        for i in link_clusters_src[d.from_segment.name]:
                            repl.append(i)
                    except(KeyError):
                        pass
                for i in repl:
                    print(str(d).replace(d.from_segment.name,'%s_%s' % (d.from_segment.name,i)))
                    try:
                        g.add_line(str(d).replace(d.from_segment.name,'%s_%s' % (d.from_segment.name,i)))
                    except(gfapy.error.NotUniqueError):
                        pass
                #print(d.from_segment.name)
                #link_clusters_src[d.from_segment.name]


def test(g):
    #gfapy.Gfa.to_file(g, gfa_transformed)
    repeat=False
    #changed = False
    #print("NEW ERA")
    #for edge in ['edge_673', 'edge_677', 'edge_806_160', 'edge_806_461', 'edge_807_113', 'edge_807_140', 'edge_807_179', 'edge_807_249']:
    for edge in g.segment_names:
        changed = clear_links2(edge)
        #print (changed)
        if changed==True:
            repeat=True

    if repeat ==True:
        test(g)


G=gfa_to_nx(g)

try:
    all_data = np.load("output/all_data.npy" , allow_pickle='TRUE').item()
except(FileNotFoundError):
    pass
for i in range(0, len(edges)):
    graph_create_unitigs(i)
    np.save("output/all_data.npy", all_data)
for i in range(0, len(edges)):
    graph_link_unitigs(i,G)
gfapy.Gfa.to_file(g,gfa_transformed)


for ed in g.segments:
    if ed.name in remove_clusters:
        g.rm(ed)
        print(ed.name)
for link in g.dovetails:
    if link.to_segment in remove_clusters or link.from_segment in remove_clusters:
        g.rm(link)

gfapy.Gfa.to_file(g,gfa_transformed)

#simplify_links(g)


gfapy.Gfa.to_file(g,gfa_transformed1)
gfapy.GraphOperations.merge_linear_paths(g)
gfapy.Gfa.to_file(g,gfa_transformed2)


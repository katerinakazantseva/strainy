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

g = gfapy.Gfa.from_file(gfa)
stats = open('output/stats_clusters.txt', 'a')
stats.write("Edge" + "\t" + "Fill Clusters" + "\t" + "Full Paths Clusters" + "\n")
stats.close()



def add_child_edge(edge, clN, g, cl, SNP_pos, data, left, righ,cons):
    f = g.segments
    seq=[]
    cl_consensuns = cluster_consensuns(cl, clN, SNP_pos, data, cons)
    for i in f:
        if i == edge:
            seq = i.sequence
            seq = list(seq)
            for key, val in cl_consensuns[clN].items():
                try:
                    seq[int(key)-1] = val
                except (ValueError, KeyError):
                    continue
    seq = ''.join(seq)
    seq = seq[left:righ]

    if len(seq)==0:
        remove_zeroes.append("S\t%s_%s\t*" % (edge, clN))
    if len(seq)>0:
        #g.add_line("S\tNEW\t*")
        g.add_line("S\t%s_%s\t*" % (edge, clN))
        f = g.segments
        for i in f:

    #if i == "NEW":
            if i == "%s_%s" % (edge, clN):
                new_line = i
                new_line.name = str(edge) + "_" + str(clN)
                new_line.sid = str(edge) + "_" + str(clN)
                new_line.sequence = seq
                new_line.dp = cons[clN]["Cov"]  # coverage
                print("edge added:" + str(new_line.name))

def remove_link(fr,fr_or, to, to_or):
    print("remove line: "+str(fr)+str(fr_or)+" to "+ str(to)+str(to_or))
    for i in g.dovetails:
        if i.from_segment==fr and i.from_orient==fr_or and i.to_segment==to and i.to_orient==to_or:
            g.rm(i)


def build_paths_graph(SNP_pos, cl, cons,full_clusters, data,ln, full_paths_roots, full_paths_leafs):
    M = build_adj_matrix_clusters(cons, SNP_pos, cl, False)
    M = change_w(M, 1)
    G = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
    G.remove_edges_from(list(nx.selfloop_edges(G)))

    for e in G.edges():
        first_cl=e[0]
        second_cl=e[1]
        intersect = set(range(cons[first_cl]["Start"], cons[first_cl]["Stop"])).intersection(
            set(range(cons[second_cl]["Start"], cons[second_cl]["Stop"])))
        G[e[0]][e[1]]['weight'] = len(intersect)
    for full_cluster in full_clusters:
        if strong_tail(full_cluster, cl, ln, "root", data)==True and strong_tail(full_cluster, cl, ln, "leaf", data)==True:
            G.remove_node(full_cluster)
            full_paths_roots.remove(full_cluster)
            full_paths_leafs.remove(full_cluster)
            #print("------------------")
    #print(G)
    path_remove = []
    G=remove_nested(G, cons)
    G = remove_not_strong_tails(G, cl, cons, data, ln, full_paths_roots, full_paths_leafs)

    for node in G.nodes():
        neighbors = nx.all_neighbors(G, node)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G, node, neighbor):
                if len(n_path) == 3:
                    path_remove.append(n_path)

    for n_path in path_remove:
        #print(n_path)
        try:
            G.remove_edge(n_path[0], n_path[2])

        except:
            continue
    return (G)

def remove_nested(G, cons):
    nodes=list(G.nodes())
    for node in nodes:
        try:
            neighbors = nx.all_neighbors(G, node)
            for neighbor in list(neighbors):
                if cons[node]["Start"]<cons[neighbor]["Start"] and cons[node]["Stop"]>cons[neighbor]["Stop"]:
                    try:
                        G.remove_node(neighbor)
                    except:
                        continue
        except:
            continue
    return (G)


def paths_graph_add_vis(edge,cons, SNP_pos, cl,full_paths_roots,full_paths_leafs):
    M = build_adj_matrix_clusters(cons, SNP_pos, cl, False)
    M = change_w(M, 1)
    G_vis = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    cl_removed = []
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    path_remove = []
    G_vis=remove_nested(G_vis, cons)
    for node in G_vis.nodes():
        neighbors = nx.all_neighbors(G_vis, node)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G_vis, node, neighbor):
                if len(n_path) == 3:
                    path_remove.append(n_path)
    #for n_path in path_remove:
        #try:
            #G_vis.remove_edge(n_path[0], n_path[1])
            #cl_removed.append(n_path[1])
            #print(n_path)
        #except:
            #continue
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
    test = nx.nx_agraph.to_agraph(G_vis)

    test = str(test)
    test = test.replace('weight', 'label')
    test = gv.AGraph(test)
    test.layout(prog="neato")
    test.draw("output/graphs/full_paths_cluster_GV_graph_%s.png" % (edge))
    G_vis.remove_node("Src")
    G_vis.remove_node("Sink")
    return(cl_removed)

def find_full_paths(G, paths_roots, paths_leafs):
    #for node in G.nodes():
       # neighbors = nx.all_neighbors(G, node)
        #for neighbor in list(neighbors):
           # n_paths = nx.algorithms.all_simple_paths(G, node, neighbor)
            #for n_path in list(n_paths):
                #if len(n_path) == 3:
                    #G.remove_edge(n_path[0], n_path[1])

    paths = []
    for root in paths_roots:
        paths_nx = nx.algorithms.all_simple_paths(G, root, paths_leafs)
        for path in list(paths_nx):
            paths.append(path)
    return (paths)


def add_link(fr, fr_or, to, to_or,w):
    link = 'L	%s	%s	%s	%s	0M	xx:Z:ex_%s' % (fr, fr_or, to, to_or, w)
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
                    w=G[path[i]][path[i+1]]['weight']
                    str='L	first_edge	+	second_edge	+	0M	xx:Z:in_%s' % w
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
    cut_l = {}
    cut_r = {}

    for path_cluster in set(path_cl):
        if path_cluster in paths_roots:
            cut_l[path_cluster] = 0
            cut_r[path_cluster] = None
        elif path_cluster in paths_leafs:
            cut_r[path_cluster] = ln - 1
            cut_l[path_cluster] = None
        else:
            cut_l[path_cluster] = None
            cut_r[path_cluster] = None
    while None in cut_l.values():
        for member in cut_l.keys():
            if cut_l[member] != None and (cut_r[member] == None or member in paths_leafs):
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
    f = g.segments
    cov=0
    len_cl=[]
    for i in othercl:
        cov=cov+cons[i]["Cov"]*(cons[i]["Stop"]-cons[i]["Start"])
        for i in range(cons[i]["Start"],cons[i]["Stop"]):
            len_cl.append(i)
    if (len(set(len_cl))/ln)<0.8 and len(clusters)-len(othercl)!=0:
        remove_clusters.append(edge)
    cov=cov/ln
    for i in f:
        if i == edge:
            i.dp =round(cov)
    return(cov)

def change_sec(g,edge, othercl, cl,SNP_pos, data):
    temp={}
    other_cl=cl
    for cluster in othercl:
        other_cl.loc[cl['Cluster']==cluster, "Cluster"] = "OTHER_%s" %edge
    cl_consensuns = cluster_consensuns(other_cl, "OTHER_%s" %edge, SNP_pos, data, temp)
    for i in g.segments:
        if i == edge:
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
            i.sequence=seq




def to_neighbours(g,edge,orient):
    to_ng=[]
    for i in g.dovetails:
        if orient=='+':
            if i.from_segment.name==edge and i.from_orient=='+':

                nei=[i.to_segment.name,i.to_orient]
                to_ng.append(nei)
        if orient == '-':
            if i.from_segment.name==edge and i.from_orient=='-':
                nei=[i.to_segment.name,i.to_orient]
                to_ng.append(nei)
    return (to_ng)



def from_neighbours(g,edge, orient):
    from_ng=[]

    for i in g.dovetails:
        if orient == '+':
            if i.to_segment.name==edge and i.to_orient=='+':
                nei = [i.from_segment.name, i.from_orient]
                from_ng.append(nei)
        if orient == '-':
            if i.to_segment.name == edge and (i.to_orient=='-'):
                nei = [i.from_segment.name, i.from_orient]
                from_ng.append(nei)
    return (from_ng)


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

def clear_links(edge):
    print('-------------------')
    print("CLEAR" + str(edge))
    print("CLEAR to + "+ str(edge))
    changed=False

    to_n=to_neighbours(g,edge,'+')
    print(to_n)
    if len(to_n)==1:
        for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
            if len(to_neighbours(g,i[0],i[1]))>1:
                remove_link(i[0],i[1], to_n[0][0],to_n[0][1])
                changed = True

    print("CLEAR to - " + str(edge))
    to_n=to_neighbours(g,edge,'-')
    print(to_n)
    if len(to_n)==1:
        #print(from_neighbours(g,to_n[0][0],to_n[0][1]))
        for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
            #print(to_neighbours(g,i[0],i[1]))
            if len(to_neighbours(g,i[0],i[1]))>1:
                remove_link(i[0],i[1], to_n[0][0],to_n[0][1])
                changed = True

    print("CLEAR from + " + str(edge))
    from_n = from_neighbours(g, edge, '+')
    print(from_n)
    if len(from_n) == 1:
        #print("len 1")
        #print(to_neighbours(g, from_n[0][0], from_n[0][1]))
        for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
            if len(from_neighbours(g, i[0], i[1])) > 1:
                remove_link(from_n[0][0], from_n[0][1],i[0], i[1])
                changed = True

    print("CLEAR from - " + str(edge))
    from_n = from_neighbours(g, edge, '-')
    print(from_n)
    if len(from_n) == 1:
        for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
            if len(from_neighbours(g, i[0], i[1])) > 1:
                remove_link(from_n[0][0], from_n[0][1],i[0], i[1])
                changed = True
    print("  ")
    return (changed)


def strong_tail(cluster, cl, ln, type, data):
    count = 0
    reads = list(cl.loc[cl['Cluster'] == cluster, 'ReadName'])
    for read in reads:
        if (data[read]["Start"] == 0 and type=="root") or (data[read]["Stop"] == ln - 1 and type=="leaf"):
            count = count + 1
    if count > 2:
        res = True
    else:
        res = False
    return (res)


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
            cons[cluster]["Start"]=2
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
            cons[cluster]["Stop"] = cons[cluster]["Stop"]-2
            #print(cons[cluster]["Stop"])
    for n_path in path_remove:
        #print(n_path)
        G.remove_edge(n_path[0], n_path[1])

    return (G)


full_cl = {}
full_paths = {}
full_paths_leafs_roots = {}
paths = {}
full_path_clusters = {}
subunits_borderline = {}
connected_subunits = {}
link_clusters = {}
link_clusters_src = {}
remove_clusters = []
remove_zeroes = []

def graph_create_unitigs(i):
    edge = edges[i]
    full_paths_roots = []
    full_paths_leafs = []
    full_clusters = []
    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
        SNP_pos = read_snp(snp, edge,bam, AF)
        data = read_bam(bam,edge,SNP_pos,clipp,min_mapping_quality,min_al_len,de_max)
        ln = int(pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4])
        clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
        cons = build_data_cons(cl, SNP_pos, data)
        for cluster in clusters:
            clStart=cons[cluster]["Start"]
            clStop = cons[cluster]["Stop"]
            if clStart == 0 and clStop == ln - 1:

                full_paths_roots.append(cluster)
                full_paths_leafs.append(cluster)
                if strong_tail(cluster, cl, ln, "root", data) == True and strong_tail(cluster, cl, ln, "leaf",
                                                                                         data) == True:
                    add_child_edge(edge, cluster, g, cl, SNP_pos, data, 0, ln - 1, cons)
                    full_clusters.append(cluster)

            elif clStart == 0:
                full_paths_roots.append(cluster)
            elif clStop == ln - 1:
                full_paths_leafs.append(cluster)
        #G = build_paths_graph(SNP_pos, cl, cons, full_clusters)
        G=build_paths_graph(SNP_pos, cl, cons, full_clusters, data, ln, full_paths_roots, full_paths_leafs)
        #G = remove_not_strong_tails(G, cl, cons,data,ln,full_paths_roots,full_paths_leafs)



        full_cl[edge] = full_clusters
        cl_removed=paths_graph_add_vis(edge,cons, SNP_pos,cl,full_paths_roots, full_paths_leafs)
        try:
            full_paths[edge] = find_full_paths(G,full_paths_roots, full_paths_leafs)
            #print(full_paths)
        except(ValueError):
            pass
        add_path_edges(edge,g,cl, data, SNP_pos, ln,full_paths, G,full_paths_roots, full_paths_leafs,cons)
        #print("paths clusters added")
        print(G)
        add_path_links(edge, full_paths[edge], G)
        #print("paths links added")
        othercl=list(set(clusters)-set(full_clusters)-set([j for i in full_paths[edge] for j in i])-set(cl_removed))
        new_cov=change_cov(g,edge,cons,ln,clusters,othercl)
        if new_cov<3 and len(clusters)-len(othercl)!=0: #PARAMETER
        #if len(clusters) - len(othercl) != 0:
            remove_clusters.append(edge)
        else:
            change_sec(g, edge, othercl, cl, SNP_pos, data)

        link_clusters[edge] = list(full_clusters) + list(
            set(full_paths_roots).intersection(set([j for i in full_paths[edge] for j in i]))) + list(
            set(full_paths_leafs).intersection(set([j for i in full_paths[edge] for j in i])))
        link_clusters_src[edge] = list(full_clusters) + list(
            set(full_paths_roots).intersection(set([j for i in full_paths[edge] for j in i])))
    except(FileNotFoundError, IndexError):
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



def graph_link_unitigs(i):
    print("CREATING NEW LINKS")
    edge = edges[i]
    #SNP_pos = read_snp(snp, edge, bam, AF)
    #data = read_bam(bam, edge, SNP_pos, clipp, min_mapping_quality, min_al_len, de_max)

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
        als = gaf.loc[gaf['ReadName'].isin(reads)]
        als = als.to_dict('split')['data']
        neighbours={}
        orient={}

        for aln in als:
            al = aln[1]
            read = aln[0]
            if re.search(">%s>.*" % edge, al):
                fr_or="+"
                to_or='+'
                n=re.sub("<.*" , "", re.sub(">.*" , "", re.sub(".*>%s>" % edge, "", al, count=0, flags=0), count=0, flags=0), count=0, flags=0)
            if re.search("<%s<.*" % edge, al):
                fr_or="-"
                to_or='-'
                n=re.sub("<.*" , "", re.sub(">.*" , "", re.sub(".*<%s<" % edge, "", al, count=0, flags=0), count=0, flags=0), count=0, flags=0)
            if re.search(">%s<.*" % edge, al):
                fr_or="+"
                to_or='-'
                n=re.sub("<.*" , "", re.sub(">.*" , "", re.sub(".*>%s<" % edge, "", al, count=0, flags=0), count=0, flags=0), count=0, flags=0)
            if re.search("<%s>.*" % edge, al):
                fr_or="-"
                to_or='+'
                n=re.sub("<.*" , "", re.sub(">.*" , "", re.sub(".*<%s>" % edge, "", al, count=0, flags=0), count=0, flags=0), count=0, flags=0)
            try:
                #common_reads[read]=[n,fr_or, to_or]
                orient[n]=[fr_or, to_or]
                neighbours[read]=n
            except(UnboundLocalError): continue

        print("neighbours")
        print(set({k for k, v in Counter(neighbours.values()).items() if v > 2}))


        #if len(set({k for k, v in Counter(neighbours.values()).items() if v > 0})) == 0:
        '''
        for read in reads:
                #print(data[read]["Lclip"])
                #print(data[read]["Rclip"])
                #print()
                for n in data[read]["Lclip"]:
                    neighbours[read]=n
                    orient[n]=['+','+']


                for n in data[read]["Lclip"]:
                    neighbours[read]=n
                    orient[n]=['-','-']'''

        for n in set({k for k, v in Counter(neighbours.values()).items() if v > 2}):
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
            n_cl = cl_n.loc[cl_n['ReadName'].isin(reads), 'Cluster']

            #print(n_cl)
            n_cl = list(
                set([x for x in list(n_cl) if Counter(list(n_cl))[x] / sum(Counter(list(n_cl)).values()) >= 0.2]))



            if len(n_cl)==0:
                try:

                    #add_link("%s_%s" % (edge, clN), fr_or, n, to_or)
                    #когда мы не находим соседей, соединяем со всеми
                    if n in remove_clusters:
                        n_cl = link_clusters_src[n]
                        #for n_cluster in link_clusters_src[n]:

                    else:
                        add_link("%s_%s" % (edge, clN), fr_or, n, to_or,w)
                    #n_cl=link_clusters_src[n]
                except (KeyError):
                    continue

            link_added=False

            for i in n_cl:
                w=Counter(list(n_cl))[i]
                try:
                    if g.try_get_segment("%s_%s" % (n, i)):
                        link_added=True
                        add_link("%s_%s" % (edge, clN), fr_or, "%s_%s" % (n, i), to_or,w)
                except(gfapy.NotFoundError):
                    continue

            if link_added==False:
                if n in remove_clusters:
                    n_cl = link_clusters_src[n]
                    # for n_cluster in link_clusters_src[n]:

                else:
                    add_link("%s_%s" % (edge, clN), fr_or, n, to_or,w)

            for i in n_cl:
                try:
                    if g.try_get_segment("%s_%s" % (n, i)):
                        link_added=True
                        add_link("%s_%s" % (edge, clN), fr_or, "%s_%s" % (n, i), to_or,w)
                except(gfapy.NotFoundError):
                    continue

def test(g):
    repeat=False
    for edge in g.segment_names:
        changed = clear_links(edge)
        if changed==True:
            repeat=True
    if repeat ==True:
        test(g)



for i in range(0, len(edges)):
    graph_create_unitigs(i)
for i in range(0, len(edges)):
    graph_link_unitigs(i)
gfapy.Gfa.to_file(g,gfa_transformed)


for ed in g.segments:
    if ed.name in remove_clusters:
        g.rm(ed)
        print("removed")
        print(ed.name)


gfapy.Gfa.to_file(g,gfa_transformed1)

test(g)
gfapy.Gfa.to_file(g,gfa_transformed2)


#gfapy.GraphOperations.merge_linear_paths(g)
#gfapy.Gfa.to_file(g,gfa_transformed)


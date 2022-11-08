import pickle
from collections import deque

import pygraphviz as gv

from cluster_postprocess import *
from build_data  import *
from params import *
from flye_consensus import FlyeConsensus

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
    try:
        seq = seq[left:righ+1]
    except TypeError:
        seq = seq[left]

    if len(seq)==0:
        remove_zeroes.append("S\t%s_%s\t*" % (edge, clN))
    if len(seq)>0:
        g.add_line("S\t%s_%s\t*" % (edge, clN))
        f = g.segments
        for i in f:

            if i == "%s_%s" % (edge, clN):
                new_line = i
                new_line.name = str(edge) + "_" + str(clN)
                new_line.sid = str(edge) + "_" + str(clN)
                new_line.sequence = seq
                new_line.dp = cons[clN]["Cov"]  # coverage
                print("edge added:" + str(new_line.name))

def remove_link(fr,fr_or, to, to_or):
    res=False
    for i in g.dovetails:
        if i.from_segment == fr and i.to_segment == to:
            g.rm(i)
            res = True
            print("remove line: " + str(i.from_segment) + str(i.from_orient) + " to " + str(i.to_segment) + str(i.to_orient))
    return (res)

def clear_links_weight(edge):
    to_lines=[]
    to_tags=[]
    fr_lines=[]
    fr_tags=[]
    for d in g.dovetails:
        if (d.from_segment == edge and d.from_orient == "+") or (d.to_segment == edge and d.from_orient == "-"):
            if d.ex != None:
                to_lines.append(d)
                to_tags.append(d.ex)
        if (d.from_segment==edge and d.from_orient == "-") or (d.to_segment == edge and d.from_orient == "+"):
            if d.ex != None:
                fr_lines.append(d)
                fr_tags.append(d.ex)
    try:
        max_index=to_tags.index(max(to_tags))
        for i in to_lines:
            if to_lines.index(i)!=max_index:
                g.rm(i)
                print("line removed: "+str(i))
    except ValueError:
        pass
    try:
        max_index = fr_tags.index(max(fr_tags))
        for i in fr_lines:
            if fr_lines.index(i)!=max_index:
                g.rm(i)
                print("line removed: "+str(i))
    except ValueError:
        pass

def build_paths_graph(SNP_pos, cl, cons,full_clusters, data,ln, full_paths_roots, full_paths_leafs, edge, flye_consensus):
    M = build_adj_matrix_clusters(cons, cl, edge, flye_consensus, False)
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
        try:
            G.remove_edge(n_path[0], n_path[2])
        except:
            continue

    return G


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
    return G


def paths_graph_add_vis(edge,cons, SNP_pos, cl,full_paths_roots,full_paths_leafs, flye_consensus):
    M = build_adj_matrix_clusters(cons, cl, edge, flye_consensus, False)
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
    return cl_removed

def find_full_paths(G, paths_roots, paths_leafs):
    paths = []
    for root in paths_roots:
        paths_nx = nx.algorithms.all_simple_paths(G, root, paths_leafs)
        for path in list(paths_nx):
            paths.append(path)
    return paths


def add_link(fr, fr_or, to, to_or,w):
    link = 'L	%s	%s	%s	%s	0M	ex:i:%s' % (fr, fr_or, to, to_or, w)
    try:
        g.add_line(link)
        print("link added from %s %s to %s %s" % (fr, fr_or, to, to_or))
    except(gfapy.NotUniqueError): pass


def add_path_links(edge, paths, G):
    for path in paths:
        for i in range(0, len(path) - 1):
                try:
                    w = G[path[i]][path[i+1]]['weight']
                    str = 'L	first_edge	+	second_edge	+	0M	ix:i:%s' % w
                    g.add_line(str.replace('first_edge', "%s_%s" % (edge, path[i])).replace(
                        'second_edge', "%s_%s" % (edge, path[i + 1])))
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
        cut_l[path_cluster] = None
        cut_r[path_cluster] = None
        if path_cluster in paths_roots:

            cut_l[path_cluster] = 0
        if path_cluster in paths_leafs:
            cut_r[path_cluster] = ln - 1

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
                                    if path[path.index(n) - 1] not in visited:
                                        R.append(path[path.index(n) - 1])
                                        Q.append(path[path.index(n) - 1])
                            except (ValueError, IndexError):
                                continue
                    else:
                        for path in paths[edge]:
                            try:
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

                L = list(set(L))
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
    return path_cl


def change_cov(g, edge, cons, ln, clusters, othercl):
    f = g.segments
    cov=0
    len_cl = []
    for i in othercl:
        cov=cov+cons[i]["Cov"]*(cons[i]["Stop"]-cons[i]["Start"])
        for i in range(cons[i]["Start"],cons[i]["Stop"]):
            len_cl.append(i)
    if (len(set(len_cl))/ln)<0.7 and len(clusters)-len(othercl)!=0:
        remove_clusters.append(edge)
    cov=cov/ln
    for i in f:
        if i == edge:
            i.dp = round(cov)
    return cov


def change_sec(g, edge, othercl, cl, flye_consensus):
    cl_copy = cl.copy()
    for cluster in othercl:
        cl_copy.loc[cl['Cluster'] == cluster, "Cluster"] = "OTHER_%s" % edge
    consensus = flye_consensus.flye_consensus("OTHER_%s" % edge, edge, cl_copy)
    g.line(edge).sequence = str(consensus['consensus'])


def to_neighbours(g,edge,orient):
    to_ng=[]
    for i in g.dovetails:
            if i.from_segment.name==edge and i.from_orient=='+':
                nei=[i.to_segment.name,i.to_orient]
                to_ng.append(nei)
            if i.to_segment.name==edge and i.to_orient=='-':
                nei=[i.from_segment.name,i.from_orient]
                to_ng.append(nei)
    return to_ng


def from_neighbours(g, edge, orient):
    from_ng=[]
    for i in g.dovetails:
        if i.to_segment.name == edge and i.to_orient=='+':
            nei = [i.from_segment.name, i.from_orient]
            from_ng.append(nei)

        if i.from_segment.name == edge and (i.from_orient=='-'):
            nei = [i.to_segment.name, i.to_orient]
            from_ng.append(nei)
    return from_ng


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
    print("CLEAR to + "+ str(edge))
    changed=False
    to_n=to_neighbours(g,edge,'+')
    print(to_n)
    if len(to_n)==1:
        print(from_neighbours(g,to_n[0][0],to_n[0][1]))
        for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
            print()
            print(i)
            print(to_neighbours(g,i[0],i[1]))
            if len(to_neighbours(g,i[0],i[1]))>1:
                remove_link(i[0],i[1], to_n[0][0],to_n[0][1])
                print("remove")
                changed = True
    print()
    print("CLEAR to - " + str(edge))
    to_n=to_neighbours(g,edge,'-')
    print(to_n)
    if len(to_n)==1:
        for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
            if len(to_neighbours(g,i[0],i[1]))>1:
                remove_link(i[0],i[1], to_n[0][0],to_n[0][1])
                changed = True

    print("CLEAR from + " + str(edge))
    from_n = from_neighbours(g, edge, '+')
    print(from_n)
    if len(from_n) == 1:
        for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
            if len(from_neighbours(g, i[0], i[1])) > 1:
                remove_link(from_n[0][0], from_n[0][1],i[0], i[1])
                changed = True
    print("CLEAR from - " + str(edge))
    from_n = from_neighbours(g, edge, '-')
    if len(from_n) == 1:
        for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
            if len(from_neighbours(g, i[0], i[1])) > 1:
                remove_link(from_n[0][0], from_n[0][1],i[0], i[1])
                changed = True
    print("  ")
    return changed


def clear_links2(edge):
    changed=False
    to_n=to_neighbours(g,edge,'+')
    if len(to_n)==1:
        for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
            if len(to_neighbours(g,i[0],i[1]))>1:
                changed=remove_link(i[0],i[1], to_n[0][0],to_n[0][1])

    from_n = from_neighbours(g, edge, '+')
    if len(from_n) == 1:
        for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
            if len(from_neighbours(g, i[0], i[1])) > 1:
                changed=remove_link(from_n[0][0], from_n[0][1],i[0], i[1])
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
    return res


def remove_not_strong_tails(G, cl,cons, data,ln, full_paths_roots, full_paths_leafs):
    path_remove=[]
    for cluster in full_paths_roots:
        if strong_tail(cluster, cl, ln, "root", data):
            neighbors = nx.all_neighbors(G, cluster)
            for neighbor in list(neighbors):
                for n_path in nx.algorithms.all_simple_paths(G, neighbor, cluster):
                    if len(n_path) == 2:
                        path_remove.append(n_path)
        else:
            full_paths_roots.remove(cluster)
            cons[cluster]["Start"]=2
    for cluster in full_paths_leafs:
        if strong_tail(cluster, cl, ln, "leaf", data):
            neighbors = nx.all_neighbors(G, cluster)
            for neighbor in list(neighbors):
                for n_path in nx.algorithms.all_simple_paths(G, cluster, neighbor):
                    if len(n_path) == 2:
                        path_remove.append(n_path)
        else:

            full_paths_leafs.remove(cluster)
            cons[cluster]["Stop"] = cons[cluster]["Stop"]-2
    for n_path in path_remove:
        try:
            G.remove_edge(n_path[0], n_path[1])
        except:
            print('tried to remove non-existant edge')

    return G


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


def graph_create_unitigs(i, flye_consensus):
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
        if len(clusters) > 1:
            for cluster in clusters:
                clStart = cons[cluster]["Start"]
                clStop = cons[cluster]["Stop"]
                if clStart == 0 and clStop == ln - 1:

                    full_paths_roots.append(cluster)
                    full_paths_leafs.append(cluster)
                    if strong_tail(cluster, cl, ln, "root", data) and strong_tail(cluster, cl, ln, "leaf", data):
                        add_child_edge(edge, cluster, g, cl, SNP_pos, data, 0, ln - 1, cons)
                        full_clusters.append(cluster)

                elif clStart == 0:
                    full_paths_roots.append(cluster)
                elif clStop == ln - 1:
                    full_paths_leafs.append(cluster)
            G = build_paths_graph(SNP_pos, cl, cons, full_clusters, data, ln, full_paths_roots, full_paths_leafs, edge,
                                flye_consensus)

            full_cl[edge] = full_clusters
            cl_removed=paths_graph_add_vis(edge,cons, SNP_pos,cl,full_paths_roots, full_paths_leafs, flye_consensus)
            try:
                full_paths[edge] = find_full_paths(G,full_paths_roots, full_paths_leafs)
            except(ValueError):
                pass
            add_path_edges(edge,g,cl, data, SNP_pos, ln,full_paths, G,full_paths_roots, full_paths_leafs,cons)
            add_path_links(edge, full_paths[edge], G)
            othercl=list(set(clusters)-set(full_clusters)-set([j for i in full_paths[edge] for j in i])-set(cl_removed))
            new_cov=change_cov(g,edge,cons,ln,clusters,othercl)
            if new_cov<6 and len(clusters)-len(othercl)!=0: #PARAMETER
                remove_clusters.append(edge)
            else:
                change_sec(g, edge, othercl, cl, flye_consensus)

            link_clusters[edge] = list(full_clusters) + list(
                set(full_paths_roots).intersection(set([j for i in full_paths[edge] for j in i]))) + list(
                set(full_paths_leafs).intersection(set([j for i in full_paths[edge] for j in i])))
            link_clusters_src[edge] = list(full_clusters) + list(
                set(full_paths_roots).intersection(set([j for i in full_paths[edge] for j in i])))
            link_clusters_sink[edge] = list(full_clusters) + list(
                set(full_paths_leafs).intersection(set([j for i in full_paths[edge] for j in i])))
        else:
            change_sec(g, edge, [clusters[0]], cl, flye_consensus)
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
    except(KeyError, UnboundLocalError):
        pass

    othercl=len(clusters)-fcN-fpN
    stats.write(edge + "\t" + str(fcN) + "\t" + str(fpN) + "\t" + str(othercl) +"\n")
    stats.close()



def graph_link_unitigs(i):
    print("CREATING NEW LINKS")
    edge = edges[i]
    clusters=[]
    try:
        clusters = link_clusters[edge]
    except KeyError:
        pass
    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
    except FileNotFoundError:
        pass
    link_unitigs=[]
    for clN in set(clusters):
        try:
            if g.try_get_segment("%s_%s" % (edge, clN)):
                link_unitigs.append(clN)
        except:
            continue

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
            n=None
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
            if n!=None:
                try:
                    orient[n]=[fr_or, to_or]
                    neighbours[read]=n
                except UnboundLocalError:
                    continue

        print("neighbours")
        print(set({k for k, v in Counter(neighbours.values()).items() if v > 1}))
        for n in set({k for k, v in Counter(neighbours.values()).items() if v > 1}):
            fr_or=orient[n][0]
            to_or=orient[n][1]
            w=1
            try:
                cl_n = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (n, I, AF), keep_default_na=False)
            except FileNotFoundError:
                add_link("%s_%s" % (edge, clN), fr_or,n, to_or,w)
                continue
            reads = []
            for k, v in neighbours.items():
                if v == n:
                    reads.append(k)
            n_cl = cl_n.loc[cl_n['ReadName'].isin(reads), 'Cluster']

            n_cl_set = list(
                set([x for x in list(n_cl) if Counter(list(n_cl))[x] / sum(Counter(list(n_cl)).values()) >= 0.2]))

            if len(n_cl_set)==0:
                try:
                    if n in remove_clusters:
                        if clN in link_clusters_sink[edge]:
                            n_cl_set = link_clusters_src[n] #
                        if clN in link_clusters_src[edge]:
                            n_cl_set = link_clusters_sink[n]

                    else:
                        add_link("%s_%s" % (edge, clN), fr_or, n, to_or,w)
                except (KeyError):
                    continue

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
                    if clN in link_clusters_sink[edge]:
                        n_cl_set = link_clusters_src[n]  #
                    if clN in link_clusters_src[edge]:
                        n_cl_set = link_clusters_sink[n]

                else:
                    add_link("%s_%s" % (edge, clN), fr_or, n, to_or,w)
            for i in n_cl_set:
                try:
                    if g.try_get_segment("%s_%s" % (n, i)):
                        add_link("%s_%s" % (edge, clN), fr_or, "%s_%s" % (n, i), to_or,w)
                except gfapy.NotFoundError:
                    continue


def test(g):
    repeat=False
    for edge in g.segment_names:
        changed = clear_links2(edge)
        if changed==True:
            repeat=True

    if repeat ==True:
        test(g)


try:
    with open(consensus_cache_path, 'rb') as f:
        print(os.getcwd())
        consensus_dict = pickle.load(f)
except:
    consensus_dict = {}

flye_consensus = FlyeConsensus(bam, gfa, 1, consensus_dict)

for i in range(0, len(edges)):
    graph_create_unitigs(i, flye_consensus)
for i in range(0, len(edges)):
    graph_link_unitigs(i)
gfapy.Gfa.to_file(g, gfa_transformed)


for ed in g.segments:
    if ed.name in remove_clusters:
        g.rm(ed)
        print(ed.name)
for link in g.dovetails:
    if link.to_segment in remove_clusters or link.from_segment in remove_clusters:
        g.rm(link)

gfapy.Gfa.to_file(g, gfa_transformed)

test(g)

gfapy.GraphOperations.merge_linear_paths(g)

flye_consensus.print_cache_statistics()
if write_consensus_cache:
    with open(consensus_cache_path, 'wb') as f:
        pickle.dump(consensus_dict, f)


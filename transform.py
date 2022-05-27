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
#gaf = pd.read_csv(gaf)


stats = open('output/stats_clusters.txt', 'a')
stats.write("Edge" + "\t" + "Fill Clusters" + "\t" + "Full Paths Clusters" + "\n")
stats.close()



def add_child_edge(edge, clN, g, cl, SNP_pos, data, left, righ,cons):
    f = g.segments
    seq=[]
    # dat = {0: "A", 1: "A", 4: "A"}
    #cons = build_data_cons(cl, SNP_pos, data)
    cl_consensuns = cluster_consensuns(cl, clN, SNP_pos, data, cons)
    #print(cl_consensuns[clN])
    for i in f:
        if i == edge:
            seq = i.sequence
            seq = list(seq)
            for key, val in cl_consensuns[clN].items():
                try:
                    seq[int(key)] = val
                except (ValueError):
                    continue

    # print(len(seq))

    seq = ''.join(seq)

    seq = seq[left:righ]


    g.add_line("S\tNEW\t*")
    f = g.segments
    for i in f:
        if i == "NEW":
            new_line = i
            new_line.name = str(edge) + "_" + str(clN)
            new_line.sid = str(edge) + "_" + str(clN)
            new_line.sequence = seq
            new_line.dp = cons[clN]["Cov"]  # coverage
            print("edge added:" + str(new_line.name))
            #print(left, righ)


def remove_link3(edge, neighbor):
    print("remove line"+str(edge)+" to "+ str(neighbor))
    edge = str(edge) + "\t"
    for i in g.edges:
        if re.search(edge, str(i)):
            if re.search(neighbor, str(i)):
                print(i)
                g.rm(i)

def remove_link(fr,fr_or, to, to_or):
    print("remove line: "+str(fr)+str(fr_or)+" to "+ str(to)+str(to_or))
    for i in g.dovetails:
        if i.from_segment==fr and i.from_orient==fr_or and i.to_segment==to and i.to_orient==to_or:
            g.rm(i)


def build_paths_graph(edge, data, SNP_pos, cl, cons):
    #cons = build_data_cons(cl, SNP_pos, data)
    M = build_adj_matrix_clusters(cons, SNP_pos, cl, False)
    M = change_w(M, 1)
    G_vis = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    for node in G_vis.nodes():
        neighbors = nx.all_neighbors(G_vis, node)
        for neighbor in list(neighbors):
            n_paths = nx.algorithms.all_simple_paths(G_vis, node, neighbor)
            for n_path in list(n_paths):
                if len(n_path) == 3:
                    print(n_path)
                    G_vis.remove_edge(n_path[0], n_path[1])

    G_vis = nx.nx_agraph.to_agraph(G_vis)
    G_vis.layout(prog="neato")
    G_vis.draw("output/graphs/full_paths_cluster_GV_graph_%s.png" % (edge))
    G = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
    return (G)

def paths_graph_add_vis(edge,G,full_paths_roots,full_paths_leafs):
    G_vis=G
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    for node in G_vis.nodes():
        neighbors = nx.all_neighbors(G_vis, node)
        for neighbor in list(neighbors):
            n_paths = nx.algorithms.all_simple_paths(G_vis, node, neighbor)
            for n_path in list(n_paths):
                if len(n_path) == 3:
                    print(n_path)
                    G_vis.remove_edge(n_path[0], n_path[1])

    G_vis.add_node("Src",style='filled',fillcolor='gray',shape='square')
    G_vis.add_node("Sink",style='filled',fillcolor='gray',shape='square')
    for i in full_paths_roots:
        G_vis.add_edge("Src", i)
    for i in full_paths_leafs:
        G_vis.add_edge(i, "Sink")
    test = nx.nx_agraph.to_agraph(G_vis)
    test.layout(prog="neato")
    test.draw("output/graphs/full_paths_cluster_GV_graph_%s.png" % (edge))
    G_vis.remove_node("Src")
    G_vis.remove_node("Sink")
    return(G_vis)

def find_full_paths(G, paths_roots, paths_leafs):
    for node in G.nodes():
        neighbors = nx.all_neighbors(G, node)
        for neighbor in list(neighbors):
            n_paths = nx.algorithms.all_simple_paths(G, node, neighbor)
            for n_path in list(n_paths):
                if len(n_path) == 3:
                    G.remove_edge(n_path[0], n_path[1])
    paths = []
    for root in paths_roots:
        paths_nx = nx.algorithms.all_simple_paths(G, root, paths_leafs)
        for path in list(paths_nx):
            paths.append(path)
    print(paths)
    return (paths)


def add_link(fr, fr_or, to, to_or):
    link = 'L	%s	%s	%s	%s	0M	RC:i:42' % (fr, fr_or, to, to_or)
    try:
        g.add_line(link)
        print("link added from %s %s to %s %s" % (fr, fr_or, to, to_or))
    except(gfapy.NotUniqueError): pass
        #print("dd")


def add_path_links(edge, paths):
    str = 'L	first_edge	+	second_edge	+	0M	RC:i:42'
    for path in paths:

        for i in range(0, len(path) - 1):
            try:
                print("path added " + str(path[i]) + " " + str(path[i + 1]))
            except:
                print("path added")
            try:
                g.add_line(str.replace('first_edge', "%s_%s" % (edge, path[i])).replace('second_edge',
                                                                                        "%s_%s" % (edge, path[i + 1])))

            except(gfapy.NotUniqueError):
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

    #print(cut_l)
    #print(cut_r)
    while None in cut_l.values():
        for member in cut_l.keys():
            if cut_l[member] != None and (cut_r[member] == None or member in paths_leafs):
                #print("NEW")
                #print(member)
                Q = deque()
                L = []
                R = []
                for path in paths[edge]:

                    try:
                        L.append(path[path.index(member) + 1])
                        #print ("add L")
                        #print(path[path.index(member) + 1])
                        Q.append(path[path.index(member) + 1])
                    except (ValueError, IndexError):
                        continue
                    #print("Q1")
                    #print(Q)

                visited = []
                while Q:
                    n = Q.pop()
                    visited.append(n)
                    if n in L:
                        for path in paths[edge]:
                            try:
                                R.append(path[path.index(n) - 1])
                                if path[path.index(n) - 1] not in visited:
                                    Q.append(path[path.index(n) - 1])


                            except (ValueError, IndexError):
                                continue
                    else:
                        for path in paths[edge]:
                            try:
                                L.append(path[path.index(n) + 1])
                                if path[path.index(n) + 1] not in visited:
                                    Q.append(path[path.index(n) + 1])

                            except (ValueError, IndexError):
                                continue
                    #print("Q2")
                    #print(Q)
                l_borders = []
                r_borders = []
                for i in L:
                    #l_borders.append(int(clusters_borders[i][0]))
                    l_borders.append(int(cons[i]["Start"]))
                for i in R:
                    #r_borders.append(int(clusters_borders[i][1]))
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
                #print("cut_l")
                #print(cut_l)
                #print("cut_r")
                #print(cut_r)

    for path_cluster in set(path_cl):
        add_child_edge(edge, path_cluster, g,  cl, SNP_pos, data, cut_l[path_cluster], cut_r[path_cluster], cons)
        #print("edge added ask")
        G.remove_node(path_cluster)
    return(path_cl)


def add_parent_subedge(edge, g, left, right):
    f = g.segments
    for i in f:
        if i == edge:
            seq = i.sequence
            seq = list(seq)
    seq = ''.join(seq)
    seq = seq[left:right]

    g.add_line("S\tNEW\t*")
    f = g.segments
    for i in f:
        if i == "NEW":
            new_line = i
            new_line.name = "%s[%s,%s]" % (edge, left, right)
            new_line.sid = "%s[%s,%s]" % (edge, left, right)
            new_line.sequence = seq
            new_line.dp = 35  # coverage
            print("edge added:" + str(new_line.name))
            print(left, right)


def change_cov(g,edge,cons,ln,othercl):
    print("coverage")
    f = g.segments
    cov=0
    for i in othercl:
        cov=cov+cons[i]["Cov"]*(cons[i]["Stop"]-cons[i]["Start"])
        print(cov)
    cov=cov/ln
    for i in f:
        if i == edge:
            i.dp =round(cov)
    print("coverage="+str(cov))
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
                    print(int(key))
                    print(val)
                    seq[int(key)] = val
                except (ValueError):
                    continue
            seq = ''.join(seq)
            i.sequence=seq
            print(seq)





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
                #print(i)
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




full_cl = {}
full_paths = {}
full_paths_leafs_roots = {}
paths = {}
full_path_clusters = {}
subunits_borderline={}
connected_subunits={}
link_clusters={}
remove_clusters=[]


def transform_graph(i):
    edge = edges[i]
    print(edge)
    full_paths_roots = []
    full_paths_leafs = []
    full_clusters = []

    full_paths_leafs_roots = {}
    # CHECK if FULL
    # cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
        #snp = None

        SNP_pos = read_snp(snp, edge,bam, AF)
        data = read_bam(bam,edge,SNP_pos,clipp,min_mapping_quality,min_al_len,de_max)
        ln = int(pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4])
        clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
        cons = build_data_cons(cl, SNP_pos, data)

        G = build_paths_graph(edge, data, SNP_pos, cl, cons)


        for cluster in clusters:

            clStart=cons[cluster]["Start"]
            clStop = cons[cluster]["Stop"]
            if clStart == 0 and clStop == ln - 1:
                full_clusters.append(cluster)
                add_child_edge(edge, cluster, g, cl, SNP_pos, data, 0, ln - 1, cons)
                #full_paths_roots.append(cluster)
                #full_paths_leafs.append(cluster)

            elif clStart == 0:
                full_paths_roots.append(cluster)

            elif clStop == ln - 1:
                full_paths_leafs.append(cluster)


        full_cl[edge] = full_clusters
        paths_graph_add_vis(edge,G, full_paths_roots, full_paths_leafs)


        try:

            full_paths[edge] = find_full_paths(G,full_paths_roots, full_paths_leafs)

            print("full paths added")
        except(ValueError):
            pass

        '''
        lis = list(nx.topological_sort(nx.line_graph(G)))

        first = []
        last = []
        parent_borders = []
        for i in lis:
            first.append(i[0])
            last.append(i[1])

        for i in lis:

            if last.count(i[0]) == 0:
                paths_roots.append(i[0])
            if first.count(i[1]) == 0:
                paths_leafs.append(i[1])



        try:
            paths[edge] = find_full_paths(G, paths_roots, paths_leafs)

        except(ValueError):
            pass
    '''
        #print(set(path_cl))
        add_path_edges(edge,g,cl, data, SNP_pos, ln,full_paths, G,full_paths_roots, full_paths_leafs,cons)
        print("paths clusters added")
        add_path_links(edge, full_paths[edge])
        print("paths links added")
        #edges_to_remove.append(edge)
        gfapy.Gfa.to_file(g,
                          "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/flye_3ecoli_sim_noalt_haplo/assembly_graph_trans.gfa")
        #full_path_clusters[edge] = set(path_cl)

        #OTHER CLUSTERS
        othercl=list(set(clusters)-set(full_clusters)-set([j for i in full_paths[edge] for j in i]))
        #print(othercl)

        new_cov=change_cov(g,edge,cons,ln,othercl)
        #print(new_cov)
        if new_cov<3: #PARAMETER
            remove_clusters.append(edge)
           # for ed in g.segments:
                #if ed.name in remove_clusters :
                    #if cut(ed) == False:
                        #g.rm(ed)
        else:
            change_sec(g, edge, othercl, cl, SNP_pos, data)
       # else:#change snips in initial edge

        link_clusters[edge] = list(full_clusters) + list(
            set(full_paths_roots).intersection(set([j for i in full_paths[edge] for j in i]))) + list(
            set(full_paths_leafs).intersection(set([j for i in full_paths[edge] for j in i])))
        #link_clusters[edge] = list(full_clusters + list(set([j for i in full_paths[edge] for j in i])))

    except(FileNotFoundError, IndexError):
        pass
        clusters = []
    # continue


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
    #edges_to_remove.append(edge)
    #if othercl==0 and k!=1:
        #edges_to_remove.append(edge)
    stats.write(edge + "\t" + str(fcN) + "\t" + str(fpN) + "\t" + str(othercl) +"\n")
    stats.close()

    # print(len(full_cl))
    # print(full_paths)
    # print(full_paths_roots,full_paths_leafs)
def transform_graph1(i):
    print(link_clusters)
    print("CREATING NEW LINKS")
    # for edge in edges:
    edge = edges[i]
    clusters=[]
    try:
        clusters = link_clusters[edge]
    except(KeyError):
        pass

    '''
    try:
        clusters = full_cl[edge]
    except(KeyError):pass

    try:
        clusters=clusters+full_paths_leafs+full_paths_roots
    except(KeyError):
        pass

    try:
        clusters = clusters + list(full_path_clusters[edge])
    except:
        pass

    try:
        clusters = clusters +list(paths_leafs)+list(paths_roots)
    except:
        pass
    # continue
    '''
    #print(clusters)


    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
    except(FileNotFoundError):
        pass

    new_links = []
    for clN in set(clusters):
        print("START")
        print("%s_%s" % (edge, clN))
        reads = list(cl.loc[cl['Cluster'] == clN, 'ReadName'])
        als = gaf.loc[gaf['ReadName'].isin(reads)]
        als = als.to_dict('split')['data']
        common_reads = {}
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

        #print(neighbours.values())

        for n in set({k for k, v in Counter(neighbours.values()).items() if v > 2}):
            fr_or=orient[n][0]
            to_or=orient[n][1]
            try:
                cl_n = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (n, I, AF), keep_default_na=False)
            except(FileNotFoundError):
                add_link("%s_%s" % (edge, clN), fr_or,n, to_or)
                continue
            reads = []
            for k, v in neighbours.items():
                if v == n:
                    reads.append(k)
            n_cl = cl_n.loc[cl_n['ReadName'].isin(reads), 'Cluster']
            n_cl = list(
                set([x for x in list(n_cl) if Counter(list(n_cl))[x] / sum(Counter(list(n_cl)).values()) >= 0.2]))
            #print(n_cl)
            #if len(full_paths[n])==0 and len(full_cl[n]==0):
                #add_link("%s_%s" % (edge, clN), fr_or, n, to_or)
            if len(n_cl)==0:
                try:
                    n_cl=link_clusters[n]
                except (KeyError):
                    continue
            #print(n_cl)
            link_added=False
            for i in n_cl:
                try:
                    g.try_get_segment("%s_%s" % (n, i))
                    link_added=True
                    add_link("%s_%s" % (edge, clN), fr_or, "%s_%s" % (n, i), to_or)
                except(gfapy.NotFoundError):
                    continue
            if link_added==False:
                add_link("%s_%s" % (edge, clN), fr_or, n, to_or)


        '''for aln in als:
            #print(aln)
            al = aln[1]
            read = aln[0]
            if re.search(".*%s<" % edge, al):
                h = re.sub("<.*", "",
                           re.sub(">.*", "", re.sub(".*%s<" % edge, "", al, count=0, flags=0), count=0, flags=0),
                           count=0,
                           flags=0)
                if len(h) > 0:
                    head[read] = h

            if re.search(".*>%s.*" % edge, al):
                h = re.sub(".*>", "",
                           re.sub(".*<", "", re.sub(">%s.*" % edge, "", al, count=0, flags=0), count=0, flags=0),
                           count=0,
                           flags=0)
                if len(h) > 0:
                    head[read] = h

            if re.search(".*%s>" % edge, al):
                t = re.sub("<.*", "",
                           re.sub(">.*", "", re.sub(".*%s>" % edge, "", al, count=0, flags=0), count=0, flags=0),
                           count=0,
                           flags=0)
                if len(t) > 0:
                    tail[read] = t
            if re.search(".*<%s.*" % edge, al):
                t = re.sub(".*>", "",
                           re.sub(".*<", "", re.sub("<%s.*" % edge, "", al, count=0, flags=0), count=0, flags=0),
                           count=0,
                           flags=0)
                if len(t) > 0:
                    tail[read] = t

            #print(head.values())
            #print(tail.values())

        for h in set({k for k, v in Counter(head.values()).items() if v > 2}):
            print(h)
            try:
                cl_h = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (h, I, AF), keep_default_na=False)
            except(FileNotFoundError):
                # print("break")
                add_link(h, "%s_%s" % (edge, clN))
                new_links.append(h)
                break
            reads = []
            for k, v in head.items():
                if v == h:
                    reads.append(k)
            h_cl = cl_h.loc[cl_h['ReadName'].isin(reads), 'Cluster']
            h_cl = list(set([x for x in list(h_cl) if Counter(list(h_cl))[x]/sum(Counter(list(h_cl)).values()) >= 0.2]))


            if len(h_cl) == 0:
                # print("no child")
                try:  # g.add_line("L\t%s\t+\t%s_%s\t-\t*" % (h, edge, clN))
                    # change_link(edge,"%s_%s" % (edge, clN), h)
                    add_link(h, "%s_%s" % (edge, clN))
                    new_links.append(h)


                except:
                    pass
                #print("addd links to " + str("%s" % (h)))
            for i in h_cl:
                if '%s_%s' % (h, i) not in g.segment_names:
                    # print("no child")
                    try:  # g.add_line("L\t%s\t+\t%s_%s\t-\t*" % (h, edge, clN))
                        add_link(h, "%s_%s" % (edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), h)
                        new_links.append(h)
                    except:
                        pass
                    #print("addd links to " + str("%s" % (h)))
                    # нет подребра, добавляем связь к головной ноде
                    pass
                else:
                    try:  # g.add_line("L\t%s_%s\t+\t%s_%s\t-\t*"  %(h, i,edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), h, i)
                        add_link("%s_%s" % (h, i), "%s_%s" % (edge, clN))
                        new_links.append(h)
                    except:
                        pass
                    #print("addd links to " + str("%s_%s" % (h, i)))

        for t in set({k for k, v in Counter(tail.values()).items() if v > 2}):
            print(t)
            try:
                cl_t = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (t, I, AF), keep_default_na=False)
            except(FileNotFoundError):
                add_link("%s_%s" % (edge, clN), t)
                new_links.append(t)
                break
            reads = []
            for k, v in head.items():
                if v == h:
                    reads.append(k)
            t_cl = cl_t.loc[cl_t['ReadName'].isin(reads), 'Cluster']
            t_cl = list(set([x for x in list(t_cl) if Counter(list(t_cl))[x] / sum(Counter(list(t_cl)).values()) >= 0.2]))

            if len(t_cl) == 0:
                #print("no child")
                try:  # g.add_line("L\t%s\t-\t%s_%s\t+\t*" % (t, edge, clN))
                    # change_link(edge, "%s_%s" % (edge, clN), t)
                    add_link("%s_%s" % (edge, clN), t)
                    new_links.append(t)
                except:
                    pass
               #print("addd links to " + str("%s" % (t)))
            for i in t_cl:
                if '%s_%s' % (h, i) not in g.segment_names:
                    # print("no child")
                    try:  # g.add_line("L\t%s\t-\t%s_%s\t+\t*" % (t, edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), t)
                        add_link("%s_%s" % (edge, clN), t)
                        new_links.append(t)
                    except:
                        pass
                    #print("addd links to " + str("%s" % (t)))

                    pass
                else:

                    try:  # g.add_line("L\t%s_%s\t-\t%s_%s\t+\t*" % (t, i, edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), t,i)
                        add_link("%s_%s" % (edge, clN), "%s_%s" % (t, i))
                        new_links.append(t)
                    except:
                        pass
                    #print("addd links to " + str("%s_%s" % (t, i)))
        #for ed in g.segments:
          #if ed.name in edges_to_remove:
               #g.rm(ed)
        #for i in new_links:
            #remove_link(edge, i)
        
        
        '''
#gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/flye_3ecoli_sim_noalt_haplo/assembly_graph_trans.gfa")


for i in range(0, len(edges)):
   transform_graph(i)

for i in range(0, len(edges)):
   transform_graph1(i)

def clear_links(edge):
    print(edge)
    changed=[]
    to_n=to_neighbours(g,edge,'+')
    if len(to_n)==1:
        for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
            if len(to_neighbours(g,i[0],i[1]))>1:
                remove_link(i[0],i[1], to_n[0][0],to_n[0][1])
                changed.append(to_n[0][0]) #+for
                changed.append(i[0]) #+to
                for k in from_neighbours(g,to_n[0][0],to_n[0][1]):
                    changed.append(k[0])
                for k in to_neighbours(g, i[0], i[1]):
                    changed.append(k[0])
    to_n=to_neighbours(g,edge,'-')
    print("too--")
    print(to_n)
    if len(to_n)==1:
        print(from_neighbours(g,to_n[0][0],to_n[0][1]))
        for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
            print(to_neighbours(g,i[0],i[1]))
            if len(to_neighbours(g,i[0],i[1]))>1:
                remove_link(i[0],i[1], to_n[0][0],to_n[0][1])
                changed.append(to_n[0][0])
                changed.append(i[0])
                for k in from_neighbours(g,to_n[0][0],to_n[0][1]):
                    changed.append(k[0])
                for k in to_neighbours(g, i[0], i[1]):
                    changed.append(k[0])
    from_n = from_neighbours(g, edge, '+')
    #print("from+")
    #print(from_n)
    if len(from_n) == 1:
        #print("len 1")
        #print(to_neighbours(g, from_n[0][0], from_n[0][1]))
        for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
            if len(from_neighbours(g, i[0], i[1])) > 1:
                remove_link(from_n[0][0], from_n[0][1],i[0], i[1])
                changed.append(from_n[0][0])
                changed.append(i[0])
                for k in to_neighbours(g,from_n[0][0],from_n[0][1]):
                    changed.append(k[0])
                for k in from_neighbours(g, i[0], i[1]):
                    changed.append(k[0])
    from_n = from_neighbours(g, edge, '-')
    if len(from_n) == 1:

        for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
            if len(from_neighbours(g, i[0], i[1])) > 1:
                remove_link(from_n[0][0], from_n[0][1],i[0], i[1])
                changed.append(from_n[0][0])
                changed.append(i[0])
                for k in to_neighbours(g,from_n[0][0],from_n[0][1]):
                    changed.append(k[0])
                for k in from_neighbours(g, i[0], i[1]):
                    changed.append(k[0])
    #print(changed)
    return (changed)

gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/flye_3ecoli_sim_noalt_haplo/assembly_graph_trans_31.gfa")




def test(edge):
    changed=clear_links(edge)
    for i in changed:
        test(i)


for edge in g.segment_names:
   test(edge)

for ed in g.segments:
 if ed.name in remove_clusters :
    #if cut(ed) == False:
    g.rm(ed)

gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/flye_3ecoli_sim_noalt_haplo/assembly_graph_trans_32.gfa")
import multiprocessing


def parallel(edges):
    pool = multiprocessing.Pool(1)
    pool.map(transform_graph, range(0, len(edges)))
    pool.close()


#if __name__ == "__main__":
    #parallel(edges)

# REMOVE LINKS


gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/flye_3ecoli_sim_noalt_haplo/assembly_graph_trans_33.gfa")



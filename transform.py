import csv
import networkx as nx
from build_adj_matrix import *
from cluster_postprocess import *
import pygraphviz as gv
import re
import gfapy
from collections import Counter, deque
from build_data  import *



edges = ['edge_1', 'edge_10', 'edge_100', 'edge_101', 'edge_102', 'edge_103', 'edge_104', 'edge_105', 'edge_106',
         'edge_107', 'edge_108', 'edge_109', 'edge_11', 'edge_110', 'edge_111', 'edge_112', 'edge_113', 'edge_114',
         'edge_115', 'edge_116', 'edge_12', 'edge_13', 'edge_14', 'edge_15', 'edge_16', 'edge_17', 'edge_18', 'edge_19',
         'edge_2', 'edge_20', 'edge_21', 'edge_22', 'edge_23', 'edge_24', 'edge_25', 'edge_26', 'edge_27', 'edge_28',
         'edge_29', 'edge_3', 'edge_30', 'edge_31', 'edge_32', 'edge_33', 'edge_34', 'edge_35', 'edge_36', 'edge_37',
         'edge_38', 'edge_39', 'edge_4', 'edge_40', 'edge_41', 'edge_42', 'edge_43', 'edge_44', 'edge_45', 'edge_46',
         'edge_47', 'edge_48', 'edge_49', 'edge_5', 'edge_50', 'edge_51', 'edge_52', 'edge_53', 'edge_54', 'edge_55',
         'edge_56', 'edge_57', 'edge_58', 'edge_59', 'edge_6', 'edge_60', 'edge_61', 'edge_62', 'edge_63', 'edge_64',
         'edge_65', 'edge_66', 'edge_67', 'edge_68', 'edge_69', 'edge_7', 'edge_70', 'edge_71', 'edge_72', 'edge_73',
         'edge_74', 'edge_75', 'edge_76', 'edge_77', 'edge_78', 'edge_79', 'edge_8', 'edge_80', 'edge_81', 'edge_82',
         'edge_83', 'edge_84', 'edge_85', 'edge_86', 'edge_87', 'edge_88', 'edge_89', 'edge_9', 'edge_90', 'edge_91',
         'edge_92', 'edge_93', 'edge_94', 'edge_95', 'edge_96', 'edge_97', 'edge_98', 'edge_99']
# edges=[ 'edge_10', 'edge_101',  'edge_103', 'edge_104', 'edge_105', 'edge_108', 'edge_11','edge_110', 'edge_111',  'edge_114', 'edge_12', 'edge_15',  'edge_2', 'edge_21', 'edge_22', 'edge_23','edge_28', 'edge_29', 'edge_31',  'edge_35', 'edge_36','edge_4','edge_42', 'edge_47', 'edge_48',  'edge_5','edge_54','edge_60',  'edge_62',  'edge_65', 'edge_71', 'edge_72',  'edge_74', 'edge_75', 'edge_79', 'edge_84',  'edge_88', 'edge_89', 'edge_90', 'edge_91', 'edge_92', 'edge_94', 'edge_96', 'edge_97', 'edge_98', 'edge_99']


# edges=['edge_96', "edge_43"]

edges = ['edge_99', 'edge_96', 'edge_98']


edges = ['edge_96']

R = 1
I = 1000
clipp = 100
min_mapping_quality = 20
min_base_quality = 0
min_al_len = 1000
de_max = 0.05
AF = 0.1

bam = "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/sim.3ecoliv2.bam"
stats = open('output/stats_clusters.txt', 'a')
stats.write("Edge" + "\t" + "Fill Clusters" + "\t" + "Full Paths Clusters" + "\n")
stats.close()



def add_child_edge(edge, clN, g, cl, SNP_pos, data, left, righ):
    f = g.segments
    seq=[]
    # dat = {0: "A", 1: "A", 4: "A"}
    cons = build_data_cons(cl, SNP_pos, data)
    cl_consensuns = cluster_consensuns(cl, clN, SNP_pos, data, cons)
    print(cl_consensuns[clN])
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
            new_line.dp = 35  # coverage
            print("edge added:" + str(new_line.name))
            print(left, righ)
    # gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph_transformed.gfa")


def remove_link(edge, neighbor):
    # print("remove "+edge+" "+ neighbor)
    edge = str(edge) + "\t"
    for i in g.edges:
        if re.search(edge, str(i)):
            if re.search(neighbor, str(i)):
                print(i)
                g.rm(i)


def build_paths_graph(edge, data, SNP_pos, cl):
    cons = build_data_cons(cl, SNP_pos, data)
    M = build_adj_matrix_clusters(cons, SNP_pos, cl)
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


def find_full_paths(G, paths_roots, paths_leafs):
    for node in G.nodes():
        neighbors = nx.all_neighbors(G, node)
        for neighbor in list(neighbors):
            n_paths = nx.algorithms.all_simple_paths(G, node, neighbor)
            for n_path in list(n_paths):
                if len(n_path) == 3:
                    print(n_path)
                    G.remove_edge(n_path[0], n_path[1])
    paths = []
    for root in paths_roots:
        paths_nx = nx.algorithms.all_simple_paths(G, root, paths_leafs)
        for path in list(paths_nx):
            paths.append(path)
    return (paths)


def add_link(first, second):
    str = 'L	first_edge	+	second_edge	+	0M	RC:i:42'
    try:
        g.add_line(str.replace('first_edge', first).replace('second_edge', second))
        print("line added " + first + " " + second)
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


def add_path_edges ( edge,g,cl, data, SNP_pos, ln, paths, G,paths_roots,paths_leafs, clusters_borders):
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
            if cut_l[member] != None and cut_r[member] == None:
                # print("NEW")
                # print(member)
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
                l_borders = []
                r_borders = []
                for i in L:
                    l_borders.append(int(clusters_borders[i][0]))
                for i in R:
                    r_borders.append(int(clusters_borders[i][1]))

                border = max(l_borders) + (min(r_borders) - max(l_borders)) // 2
                for i in L:
                    cut_l[i] = border

                for i in R:
                    cut_r[i] = border


    for path_cluster in set(path_cl):
        add_child_edge(edge, path_cluster, g,  cl, SNP_pos, data, cut_l[path_cluster], cut_r[path_cluster])
        print("edge added ask")
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


def add_parent_subunits(edge, g, parent_borders, ln):
    parent_subunits = []
    parent_borders.append(0)
    parent_borders.append(ln - 1)
    parent_borders=list(set(parent_borders))
    subunits_borderline[edge]= [[], []]



    for i in range(0, len(parent_borders)-1):
        add_parent_subedge(edge, g, sorted(parent_borders)[i], sorted(parent_borders)[i + 1])
        parent_subunits.append("%s[%s,%s]" % (edge, sorted(parent_borders)[i], sorted(parent_borders)[i + 1]))

        if parent_borders[i]==0:
            subunits_borderline[edge][0].append("%s[%s,%s]" % (edge, sorted(parent_borders)[i], sorted(parent_borders)[i + 1]))
        if parent_borders[i+1]==ln-1:
            subunits_borderline[edge][1].append("%s[%s,%s]" % (edge, sorted(parent_borders)[i], sorted(parent_borders)[i + 1]))
    return (parent_subunits)


g = gfapy.Gfa.from_file(
    "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph.gfa")
I = 1000
AF = 0.1
full_cl = {}
full_paths = {}
paths = {}
full_path_clusters = {}
subunits_borderline={}
connected_subunits={}


def transform_graph(i):
    edge = edges[i]
    print(edge)
    full_paths_roots = []
    full_paths_leafs = []
    full_clusters = []
    paths_roots = []
    paths_leafs = []
    part_clusters = []
    edges_to_remove=[]
    # CHECK if FULL
    # cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
        snp = None
        SNP_pos = read_snp(snp, edge,bam, AF)
        data = read_bam(bam,edge,SNP_pos,clipp,min_mapping_quality,min_al_len,de_max)
        ln = int(pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4])
        clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))

        clusters_borders = {}
        G = build_paths_graph(edge, data, SNP_pos, cl)

        for cluster in clusters:
            clStart = int(ln)
            clStop = 0
            for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
                start = int(data[read]["Start"])
                stop = int(data[read]["Stop"])
                if start < clStart:
                    clStart = start
                if stop > clStop:
                    clStop = stop

            if clStart == 0 and clStop == ln - 1:

                full_clusters.append(cluster)
                add_child_edge(edge, cluster, g, cl, SNP_pos, data, 0, ln - 1)
                G.remove_node(cluster)

            elif clStart == 0:
                full_paths_roots.append(cluster)

            elif clStop == ln - 1:
                full_paths_leafs.append(cluster)

            clusters_borders[cluster] = [clStart, clStop]
        full_cl[edge] = full_clusters
        # CHECK if SEMI-FULL

        try:

            full_paths[edge] = find_full_paths(G,full_paths_roots, full_paths_leafs)
        except(ValueError):
            pass
        # continue
        path_cl = add_path_edges(edge,g,cl, data, SNP_pos, ln,full_paths, G,full_paths_roots, full_paths_leafs,clusters_borders)
        add_path_links(edge, full_paths[edge])

        gfapy.Gfa.to_file(g,
                          "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph_transformed.gfa")
        full_path_clusters[edge] = set(path_cl)
        print(len(set(path_cl)))

    except(FileNotFoundError):
        pass
    # continue


    stats = open('output/stats_clusters.txt', 'a')
    fcN = 0
    fpN = 0


    try:
        fcN = len(full_cl[edge])
    except(KeyError):
        pass
    try:
        fpN = len(set(path_cl))
    except(KeyError):
        pass

    othercl=len(clusters)-fcN-fpN

    if othercl==0:
        edges_to_remove.append(edge)
    stats.write(edge + "\t" + str(fcN) + "\t" + str(fpN) + "\t" + str(othercl) +"\n")
    stats.close()

    # print(full_cl)
    # print(len(full_cl))
    # print(full_paths)
    # print(full_paths_roots,full_paths_leafs)

    print("CREATING NEW LINKS")
    # for edge in edges:
    clusters = []
    print()
    try:
        clusters = full_cl[edge]
    except(KeyError):pass

    try:
        clusters = clusters + list(full_path_clusters[edge])
    except:
        pass

    try:
        clusters = clusters +list(paths_leafs)+list(paths_roots)
    except:
        pass
    # continue
    print(clusters)

    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
    except(FileNotFoundError):
        pass
    # gaf = pd.read_csv("~/MT/gaf3noalt/gaf_%s.csv" % (edge))
    gaf = pd.read_csv("~/MT/gaf3noalt/gaf_common.csv")
    new_links = []
    for clN in clusters:
        print("%s_%s" % (edge, clN))
        reads = list(cl.loc[cl['Cluster'] == clN, 'ReadName'])
        als = gaf.loc[gaf['ReadName'].isin(reads)]
        als = als.to_dict('split')['data']
        head = {}
        tail = {}
        for aln in als:
            al = aln[2]
            read = aln[1]
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

        for h in set({k for k, v in Counter(head.values()).items() if v > 2}):
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
            h_cl = list(set(h_cl))
            print(h_cl)

            if len(h_cl) == 0:
                # print("no child")
                try:  # g.add_line("L\t%s\t+\t%s_%s\t-\t*" % (h, edge, clN))
                    # change_link(edge,"%s_%s" % (edge, clN), h)
                    add_link(h, "%s_%s" % (edge, clN))
                    new_links.append(h)


                except:
                    pass
                print("addd links to " + str("%s" % (h)))
            for i in h_cl:
                if '%s_%s' % (h, i) not in g.segment_names:
                    # print("no child")
                    try:  # g.add_line("L\t%s\t+\t%s_%s\t-\t*" % (h, edge, clN))
                        add_link(h, "%s_%s" % (edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), h)
                        new_links.append(h)
                    except:
                        pass
                    print("addd links to " + str("%s" % (h)))
                    # нет подребра, добавляем связь к головной ноде
                    pass
                else:
                    try:  # g.add_line("L\t%s_%s\t+\t%s_%s\t-\t*"  %(h, i,edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), h, i)
                        add_link("%s_%s" % (h, i), "%s_%s" % (edge, clN))
                        new_links.append(h)
                    except:
                        pass
                    print("addd links to " + str("%s_%s" % (h, i)))

        for t in set({k for k, v in Counter(tail.values()).items() if v > 2}):
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
            t_cl = list(set(t_cl))

            if len(t_cl) == 0:
                print("no child")
                try:  # g.add_line("L\t%s\t-\t%s_%s\t+\t*" % (t, edge, clN))
                    # change_link(edge, "%s_%s" % (edge, clN), t)
                    add_link("%s_%s" % (edge, clN), t)
                    new_links.append(t)
                except:
                    pass
                print("addd links to " + str("%s" % (t)))
            for i in t_cl:
                if '%s_%s' % (h, i) not in g.segment_names:
                    # print("no child")
                    try:  # g.add_line("L\t%s\t-\t%s_%s\t+\t*" % (t, edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), t)
                        add_link("%s_%s" % (edge, clN), t)
                        new_links.append(t)
                    except:
                        pass
                    print("addd links to " + str("%s" % (t)))

                    pass
                else:

                    try:  # g.add_line("L\t%s_%s\t-\t%s_%s\t+\t*" % (t, i, edge, clN))
                        # change_link(edge, "%s_%s" % (edge, clN), t,i)
                        add_link("%s_%s" % (edge, clN), "%s_%s" % (t, i))
                        new_links.append(t)
                    except:
                        pass
                    print("addd links to " + str("%s_%s" % (t, i)))
        for ed in g.segments:
           if ed.name in edges_to_remove:
               g.rm(ed)
        #for i in new_links:
            #remove_link(edge, i)
        gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph_transformed.gfa")


for i in range(0, len(edges)):
   transform_graph(i)



import multiprocessing


def parallel(edges):
    pool = multiprocessing.Pool(1)
    pool.map(transform_graph, range(0, len(edges)))
    pool.close()


#if __name__ == "__main__":
    #parallel(edges)

# REMOVE LINKS


gfapy.Gfa.to_file(g,
                  "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph_transformed.gfa")



import networkx as nx
import logging

from metaphase.clustering.community_detection import find_communities
from metaphase.clustering.build_adj_matrix import *
from metaphase.clustering.build_data import *


logger = logging.getLogger()


def split_cluster(cl,cluster, data,clSNP, bam, edge, R, I,only_with_common_snip=True):
    logging.info("Split cluster: " + str(cluster)+ " "+ str(only_with_common_snip))
    child_clusters = []
    reads=sorted(set(cl.loc[cl['Cluster'] == cluster,'ReadName'].values))
    if cluster==unclustered_group_N or cluster==unclustered_group_N2  or only_with_common_snip==False: #NA cluster
    #if cluster == unclustered_group_N and len(clSNP)==0:  # NA cluster
        #print("Build na matrix")
        m = build_adj_matrix(cl[cl['Cluster'] == cluster], data, clSNP, I, bam, edge, R, only_with_common_snip=False)
    else:
        m=build_adj_matrix(cl[cl['Cluster'] == cluster], data, clSNP, I, bam,edge,R)
    m = remove_edges(m, 1)
    m.columns=range(0,len(cl[cl['Cluster'] == cluster]['ReadName']))
    m.index=range(0,len(cl[cl['Cluster'] == cluster]['ReadName']))

    m = change_w(m, R)
    G_sub = nx.from_pandas_adjacency(m)
    cl_exist = sorted(set(cl.loc[cl['Cluster'] != 'NA','Cluster'].values))
    cluster_membership = find_communities(G_sub)
    clN=0
    uncl=0
    reads = cl[cl['Cluster'] == cluster]['ReadName'].values

    new_cl_id_na = cluster + split_id

    while new_cl_id_na in cl_exist:
        new_cl_id_na = new_cl_id_na + 1
    print(cluster_membership)
    print(new_cl_id_na)

    if len(set(cluster_membership.values()))>0:
        for value in set(cluster_membership.values()):
            group = [k for k, v in cluster_membership.items() if v == value]
            if len(group) > min_cluster_size:
                clN = clN + 1
                new_cl_id=new_cl_id_na+clN
                while new_cl_id in cl_exist:
                    new_cl_id=new_cl_id+1
                    child_clusters.append(new_cl_id)
                cl_exist.append(new_cl_id)
                for i in group:
                     mask = cl['ReadName'] == str(reads[i])
                     cl.loc[mask, "Cluster"] = new_cl_id

            else:
                uncl = uncl + 1
                for i in group:
                    mask = cl['ReadName'] == str(reads[i])
                    #cl.loc[mask, "Cluster"] = 'NA'

                if only_with_common_snip == True or cluster==1000000: #change it for parameter process NA or not
                    cl.loc[mask, "Cluster"] = new_cl_id_na
                    child_clusters.append(new_cl_id_na)
                else:
                    cl.loc[mask, "Cluster"] = 'NA'
                    #print('NA')'''
    if cluster==236:
        reads = cl[cl['Cluster'] == cluster]['ReadName'].values
        import numpy as np
        import matplotlib as mt
        #G_vis = nx.nx_agraph.to_agraph(G_sub)
        #G_vis.layout(prog="neato")
        #G_vis.draw("%s/graphs/test.png" % MetaPhaseArgs.output, node_size=10)
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap('viridis')
        #cmap = cmap(np.linspace(0, 1, 3))
        colors = {}
        strains=['Badread_Escherichia_coli_B1109', 'Badread_E.coli.b2207', 'Badread_Escherichia_coli_B3008']
        i=0
        cl.loc[cl['Cluster'] == 'NA', 'Cluster'] = 0
        clusters = sorted(set(cl['Cluster'].astype(int)))
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


        #for strain in strains:
            #colors[strain] = mt.colors.to_hex(cmap[i])
            #i = i + 1

        for index in cl.index:
            #cl.loc[index, 'Color'] = colors[str(cl.loc[index,'ReadName']).split(',')[0]]
            cl.loc[index, 'Color'] = colors[int(cl.loc[index, 'Cluster'])]

        cl = cl[cl['Cluster'] == 236]
        print(cl['Color'])
        print(G_sub.nodes())
        nx.draw(G_sub, nodelist=G_sub.nodes(),with_labels=True, width=0.03, node_size=10, font_size=5,node_color=cl['Color'])
        #plt.show()
        plt.suptitle(str(edge))
        plt.savefig("test.png" , format="PNG", dpi=300)
        plt.close()
    else:
        print("not splitted")

    print("cluster splitted for "+str(clN))
    return ([new_cl_id_na, clN])
"""
def build_adj_matrix_clusters_atab (cons, cl, edge, flye_consensus, only_with_common_snip=True):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA','Cluster'].values))
    sort=[]
    for k,v in cons.items():
        sort.append([k, int(v["Start"]), int(v["Stop"])])
    sorted_by_pos=[]
    for i in sorted(sort, key=lambda sort: [sort[1], sort[2]]):
        sorted_by_pos.append(i[0])
    clusters = sorted(set(sorted_by_pos) & set(clusters), key=sorted_by_pos.index)
    m = pd.DataFrame(-1, index=clusters, columns=clusters)
    for i in range(0, m.shape[1]):
        first_cl = m.index[i]
        for k in range(i+1,m.shape[1]):
            second_cl = m.index[k]
            if m[second_cl][first_cl] == -1:
                m[second_cl][first_cl] = flye_consensus.cluster_distance_via_alignment(first_cl, second_cl, cl, edge)
    return m
"""


def build_adj_matrix_clusters (edge,cons,cl,flye_consensus, only_with_common_snip=True):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA','Cluster'].values))
    try:
        clusters.remove(0)
    except:
        pass
    Y=[]
    X=[]
    Z=[]
    sort=[]

    for k,v in cons.items():
        X.append(k)
        Y.append(int(v["Start"]))
        Z.append(int(v["Stop"]))
        sort.append([k,int(v["Start"]),int(v["Stop"])])
    sorted_by_pos=[]

    for i in sorted(sort, key=lambda sort: [sort[2], sort[1]]):
        sorted_by_pos.append(i[0])
    clusters=sorted(set(sorted_by_pos) & set(clusters), key=sorted_by_pos.index)
    #print(clusters)
    m = pd.DataFrame(-1, index=clusters, columns=clusters)

    for i in range(0, m.shape[1]):
        first_cl=m.index[i]
        for k in range(i + 1, m.shape[1]):
            second_cl=m.index[k]

            if m[second_cl][first_cl] == -1:
                m[second_cl][first_cl] = distance_clusters(edge, first_cl, second_cl, cons, cl,flye_consensus, only_with_common_snip)
    return m


def join_clusters(cons, cl, R, edge, consensus, only_with_common_snip=True):
#def join_clusters(cons, cl, R, edge, only_with_common_snip=True):
    if only_with_common_snip==False:
        #M = build_adj_matrix_clusters(cons, cl, edge, consensus, False)
        M = build_adj_matrix_clusters(edge,cons, cl,consensus, False)
    else:
        #M = build_adj_matrix_clusters(cons, cl, edge, consensus)
        M = build_adj_matrix_clusters(edge,cons, cl,consensus)

    M=change_w(M,R)
    G_vis = nx.from_pandas_adjacency(M, create_using=nx.DiGraph)
    G_vis.remove_edges_from(list(nx.selfloop_edges(G_vis)))
    to_remove = []
    G_vis_before = nx.nx_agraph.to_agraph(G_vis)
    G_vis_before.layout(prog="neato")
    G_vis_before.draw("%s/graphs/cluster_GV_graph_before_remove_%s.png" % (MetaPhaseArgs.output, edge))


    path_remove=[]
    for node in G_vis.nodes():
        neighbors = nx.all_neighbors(G_vis, node)
        for neighbor in list(neighbors):
            for n_path in nx.algorithms.all_simple_paths(G_vis, node, neighbor, cutoff=5):
                if len(n_path) == 3:
                    path_remove.append(n_path)

    for n_path in path_remove:
        try:
            G_vis.remove_edge(n_path[0], n_path[2])
        except:
            continue

    lis = list(nx.topological_sort(nx.line_graph(G_vis)))
    first = []
    last = []

    for i in lis:
        first.append(i[0])
        last.append(i[1])

    for i in lis:
        if first.count(i[0]) > 1 or last.count(i[1]) > 1:
            to_remove.append(i)
    G_vis.remove_edges_from(ebunch=to_remove)
    G_vis = nx.nx_agraph.to_agraph(G_vis)
    G_vis.layout(prog="neato")
    G_vis.draw("%s/graphs/cluster_GV_graph_%s.png" % (MetaPhaseArgs.output, edge))
    G = nx.from_pandas_adjacency(M)
    for n_path in path_remove:
        try:
            G.remove_edge(n_path[0], n_path[2])
        except :
            continue
    G.remove_edges_from(ebunch=to_remove)


    nodes = list(G.nodes())
    for node in nodes:
        try:
            neighbors = nx.all_neighbors(G, node)
            for neighbor in list(neighbors):
                if cons[node]["Start"] < cons[neighbor]["Start"] and cons[node]["Stop"] > cons[neighbor]["Stop"]:
                    #G.add_edge(str(node),str(neighbor))
                    try:
                        G.remove_edge(node, neighbor)
                        logger.debug("REMOVE NESTED" + str(neighbor))
                    except:
                        continue
        except:
            continue

    groups = list(nx.connected_components(G))

    for group in groups:
        if len(group) > 1:
            for i in range(0, len(group)):
                cl.loc[cl['Cluster'] == list(group)[i], 'Cluster'] = int(list(group)[0])+10000000
    return cl

def split_all(cl, cluster, data, cons,bam, edge, R, I, SNP_pos,reference_seq):
    #i=i+1
    #cl.to_csv("%s/clusters/%s.csv" % (MetaPhaseArgs.output,i))
    if cons[cluster]["Strange"] == 1:
        clSNP = cons[cluster]["clSNP"]
        res = split_cluster(cl, cluster, data, clSNP, bam, edge, R, I)
        new_cl_id_na=res[0]
        clN =res[1]
        cluster_consensuns(cl, new_cl_id_na, SNP_pos, data, cons, edge, reference_seq)

        if clN!=0: #if clN==0 we dont need split NA cluster
            split_cluster(cl, new_cl_id_na, data, cons[new_cl_id_na]["clSNP"], bam, edge, R, I, False)
        clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA', 'Cluster'].values))

        if clN==1: #STOP LOOP IF EXIST
            logging.info("Probably cluster has not been splitted : " + str(new_cl_id_na + clN))
            cluster_consensuns(cl, new_cl_id_na + clN, SNP_pos, data, cons, edge, reference_seq)

        print(clusters)
        for cluster in clusters:
            if cluster not in cons:
                cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge, reference_seq)
                split_all(cl, cluster, data, cons,bam, edge, R, I, SNP_pos,reference_seq)

def split_all2(cl, cluster, data, cons,bam, edge, R, I, SNP_pos,reference_seq):
    if cons[cluster]["Strange2"] == 1:
        clSNP = cons[cluster]["clSNP2"]
        res = split_cluster(cl, cluster, data, clSNP, bam, edge, R, I)
        new_cl_id_na=res[0]
        clN =res[1]
        cluster_consensuns(cl, new_cl_id_na, SNP_pos, data, cons, edge, reference_seq)

        if clN!=0: #if clN==0 we dont need split NA cluster
            split_cluster(cl, new_cl_id_na, data, cons[new_cl_id_na]["clSNP2"], bam, edge, R, I, False)
        clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA', 'Cluster'].values))

        if clN==1: #STOP LOOP IF EXIST
            logging.info("Probably cluster has not been splitted : " + str(new_cl_id_na + clN))
            cluster_consensuns(cl, new_cl_id_na + clN, SNP_pos, data, cons, edge, reference_seq)

        print(clusters)
        for cluster in clusters:
            if cluster not in cons:
                cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge, reference_seq)
                split_all(cl, cluster, data, cons,bam, edge, R, I, SNP_pos,reference_seq)



def postprocess(bam, cl, SNP_pos, data, edge, R, I, flye_consensus):
    reference_seq = read_fasta_seq(MetaPhaseArgs.fa, edge)

    cons = build_data_cons(cl, SNP_pos, data, edge, reference_seq)
    cl.to_csv("%s/clusters/1.csv" % MetaPhaseArgs.output)
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA','Cluster'].values))
    cl.loc[cl['Cluster'] == 1000000, 'Cluster'] = 'NA'
    #clSNP = cons[1000000]["clSNP"]
    #split_cluster(cl, 1000000, data, clSNP, bam, edge, R, I, False)
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA', 'Cluster'].values))


    for cluster in clusters:
        split_all(cl, cluster, data, cons,bam, edge, R, I, SNP_pos,reference_seq)
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA', 'Cluster'].values))
    #cl.to_csv("%s/clusters/2.csv" % MetaPhaseArgs.output)

    logging.info("Split2")
    for cluster in clusters:
        split_all2(cl, cluster, data, cons,bam, edge, R, I, SNP_pos,reference_seq)
    cl.to_csv("%s/clusters/3.csv" % MetaPhaseArgs.output)


    cl.loc[cl['Cluster'] == 'NA', 'Cluster'] = 1000000
    cluster_consensuns(cl, 1000000, SNP_pos, data, cons, edge, reference_seq)
    clSNP = cons[1000000]["clSNP"]
    #split_cluster(cl, 1000000, data, clSNP, bam, edge, R, I, True)
    split_all(cl, 1000000, data, cons, bam, edge, R, I, SNP_pos, reference_seq)
    cl = cl[cl['Cluster'] != 'NA']
    cl = cl[cl['Cluster'] != 1000000]
    cl.to_csv("%s/clusters/4.csv" % MetaPhaseArgs.output)
    counts = cl['Cluster'].value_counts(dropna=False)
    cl = cl[~cl['Cluster'].isin(counts[counts < 6].index)]  # change for cov*01.

    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA', 'Cluster'].values))
    for cluster in clusters:
        if cluster not in cons:
            cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge, reference_seq)

    cl = join_clusters(cons, cl, R, edge, flye_consensus)
    cl.to_csv("%s/clusters/5.csv" % MetaPhaseArgs.output)
    #cons = build_data_cons(cl, SNP_pos, data, edge, reference_seq)
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA', 'Cluster'].values))
    for cluster in clusters:
        if cluster not in cons:
            cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge, reference_seq)
    print(sorted(set(cl.loc[cl['Cluster'] != 'NA', 'Cluster'].values)))
    cl = join_clusters(cons, cl, R, edge, flye_consensus, False)
    return(cl)



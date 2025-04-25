import pysam
import networkx as nx
import numpy as np
from strainy.color_bam import color
import multiprocessing
import gfapy
import matplotlib.pyplot as plt
import matplotlib as mt
import logging
import pandas as pd
from strainy.clustering.community_detection import find_communities
from strainy.clustering import build_adj_matrix as matrix
from strainy.clustering.cluster_postprocess2 import pairs
from strainy.clustering.cluster_postprocess2 import pairs_parallel
from strainy.clustering import build_data as build_data
from strainy.params import StRainyArgs, init_global_args_storage
from strainy.params import *
import networkit as nk
from strainy.reversed import minimize
import strainy.graph_operations.gfa_ops as gfa_ops
from strainy.flye_consensus import FlyeConsensus
from strainy.overlap import StrainyOverlap2
from multiprocessing import Pool, cpu_count
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
import os, gzip, subprocess
from Bio import SeqIO
from Bio import SeqIO
import subprocess
import gzip
import os

logger = logging.getLogger()

def all_main(args):
    #logging.info(StRainyArgs)
    logger.info("Starting phasing")
    dirs = ("%s/vcf/" % StRainyArgs().output_intermediate,
            "%s/clusters/" % StRainyArgs().output_intermediate,
            "%s/bam/merged" % StRainyArgs().output_intermediate,
            "%s/bam/clusters" % StRainyArgs().output_intermediate,
            "%s/flye_inputs" % StRainyArgs().output_intermediate,
            "%s/graphs" % StRainyArgs().output_intermediate,
            "%s/flye_outputs/asm" % StRainyArgs().output_intermediate
)
    debug_dirs = ("%s/graphs/" % StRainyArgs().output_intermediate,
                  "%s/adj_M/" % StRainyArgs().output_intermediate
    )

    for dir in dirs:
        os.makedirs(dir, exist_ok=True)
    logger.info("### Reading SNPs...")
    snp_pos = build_data.read_snp2(StRainyArgs().snp, StRainyArgs().bam, StRainyArgs().AF)
    input_graph = gfapy.Gfa.from_file(args.gfa_ref)
    edges = input_graph.segment_names
    #edges = ["edge_192","edge_188"]
    #edges=["edge_457","edge_466","edge_467","edge_468","edge_469","edge_470","edge_471","edge_472","edge_473"]
    #print(edges)
    unphased_edges = []
    for edge in edges:
        try:
            snp_pos[edge]
        except KeyError:
            snp_pos[edge] = []
            unphased_edges.append(edge)
    #print(snp_pos)
    logger.info("### Reading Reads...")
    #data = build_data.read_bam2(StRainyArgs().bam, edges, snp_pos, min_mapping_quality, min_base_quality, min_al_len,
                                #de_max[StRainyArgs().mode])

    data = build_data.read_bam2_parallel(StRainyArgs().bam, edges, snp_pos, min_mapping_quality, min_base_quality, min_al_len,
                                de_max[StRainyArgs().mode],n_processes=4)
    cl = pd.DataFrame(columns=['ReadName', 'Cluster', 'Coordinates'])

    for key, value in data.items():
        coord = {k: [v["Start"], v["End"]] for k, v in value.items()}
        row = pd.DataFrame({'ReadName': [key], 'Cluster': ['NA'], 'Coordinates': [coord]})
        cl = pd.concat([cl, row])
    cl = cl.reset_index(drop=True)
    logger.info("### Creating connection graph...")
    g="nx" #"nk"
    if g=="nx":
        m = matrix.build_adj_matrix2(edges, data, snp_pos, I, StRainyArgs().bam, 0)
        logger.info("matrix done")
        m = matrix.remove_edges(m, 0)
        m.columns = range(0, len(list(data.keys())))
        m.index = range(0, len(list(data.keys())))
        G = gfa_ops.from_pandas_adjacency_notinplace(matrix.change_w(m.transpose(), 0))
        """
        #TODO add NA treathment
        reads = sorted(set(cl.loc[cl["Cluster"] == "NA", "ReadName"].values))
        data_filtered = {key: data[key] for key in reads}
        m = matrix.build_adj_matrix2(edges, data_filtered, snp_pos, I, StRainyArgs().bam, 0, False)
        m = matrix.remove_edges(m, 0)
        m.columns = range(0, len(list(data_filtered.keys())))
        m.index = range(0, len(list(data_filtered.keys())))
        G = gfa_ops.from_pandas_adjacency_notinplace(matrix.change_w(m.transpose(), 0))
        indexes = cl[cl['ReadName'].isin(reads)].index.values.tolist()
        """
        #parrallel test
        """
        G = matrix.build_graph_from_data(
            edges, data, snp_pos, I, StRainyArgs().bam, 0,
            only_with_common_snip=True,
            use_memmap=True,
            n_processes=2,
            weight_threshold=0.0)
        """
        logger.info("### Searching clusters...")
        cluster_membership = find_communities(G)

        clN = 0
        uncl = 0
        cl_exist = []
        for value in set(cluster_membership.values()):
            group = [k for k, v in cluster_membership.items() if v == value]
            if len(group) > 3:
                clN = clN + 1
                cl.loc[group, 'Cluster'] = value
                cl_exist.append(value)
            else:
                uncl = uncl + 1
        print("UNCLUSTERED:" + str(uncl))

    #TODO check clustering
    if g == "nk":
        G1=matrix.build_graph_from_data_networKit(data, snp_pos, edges, 0, only_with_common_snip=True,
                                    weight_threshold=0, n_processes=None, chunk_size=10000)
        logger.info("### Searching clusters...")
        cluster_membership1 = nk.community.detectCommunities(G1, algo=nk.community.PLP(G1))
        clN = 0
        uncl = 0
        cl_exist = []
        for value in cluster_membership1.getSubsetIds():
            group = list(cluster_membership1.getMembers(value))
            if len(group) > 3:
                clN = clN + 1
                cl.loc[group, 'Cluster'] = value
                cl_exist.append(value)
            else:
                uncl = uncl + 1
        print("UNCLUSTERED:" + str(uncl))


    cl.to_csv(
        "%s/clusters/clusters_before_splitting_%s_%s.csv" % (StRainyArgs().output_intermediate, I, StRainyArgs().AF))
    color_bam2(edges)




    empty_consensus_dict = {}
    default_manager = multiprocessing.Manager()
    shared_flye_consensus = FlyeConsensus(StRainyArgs().bam, args.fasta_ref, 4,
                                          empty_consensus_dict, default_manager)

    # TODO add code for reversed edges
    reversed_edges=minimize(input_graph)["neg"]
    print(reversed_edges)


    logger.info("### Creating overlap graph...")
    ov_nodes = sorted(set(cl.loc[cl["Cluster"] != "NA", "Cluster"].values))
    ov_edges = pairs(edges, cl, shared_flye_consensus, data, snp_pos, reversed_edges, only_with_common_snip=True)
    """
    #ov_edges = pairs_parallel(
        #edges, cl, shared_flye_consensus, data, snp_pos, reversed_edges,
        #output_intermediate=StRainyArgs().output_intermediate,
        #only_with_common_snip=True)
    """
    over2 = StrainyOverlap2(ov_nodes, ov_edges)
    over2.remove_transitive()
    merged = over2.merge_unbrunching()
    for pair in merged:
        cl.loc[cl["Cluster"] == int(list(pair)[1]), "Cluster"] = int(list(pair)[0])
    over2.vis()
    logger.info("### Creating gfa graph...")

    args = StRainyArgs()
    run_clusters_asm_parallel(
        ov_nodes,
        cl,
        fq_path=args.fq,
        output_dir=args.output_intermediate,
        flye_path="/Users/ekaterina.kazantseva/strainy2/strainy/submodules/Flye/bin/flye",
        n_processes=4
    )
    ov_nodes = over2.nodes
    ov_edges = over2.edges
    graph = gfapy.Gfa()
    for i in ov_nodes:
        #clusters_asm(i, cl)
        for record in SeqIO.parse(f"{StRainyArgs().output_intermediate}/flye_outputs/asm/{i}/asm/assembly.fasta",
                                  "fasta"):
            seq = record.seq
        gfa_ops.add_edge(graph, str(i), 1, seq)
    for i in ov_edges:
        gfa_ops.add_link(graph, str(i[0]), "+", str(i[1]), "+", 1)
    gfapy.Gfa.to_file(graph, "%s/final.gfa" % (StRainyArgs().output))
    color_bam2(edges)
    logger.info("Done")





def clusters_asm_par(cluster, cl, fq_path, output_dir, flye_path):
    flye_dir = f"{output_dir}/flye_outputs/asm/{cluster}/asm"
    cluster_fastq = f"{output_dir}/flye_outputs/asm/{cluster}/cluster_reads.fastq"

    os.makedirs(os.path.dirname(cluster_fastq), exist_ok=True)

    out_file = open(cluster_fastq, "w")
    extracted_reads = cl.loc[cl["Cluster"] == cluster]["ReadName"].to_numpy()

    with gzip.open(fq_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if record.id in extracted_reads:
                SeqIO.write(record, out_file, "fastq")
    out_file.close()

    flye_cmd = f"{flye_path} --pacbio-hifi {cluster_fastq} -o {flye_dir}"
    subprocess.run(flye_cmd, shell=True, stderr=subprocess.DEVNULL)



def run_clusters_asm_parallel(ov_nodes, cl, fq_path, output_dir, flye_path, n_processes=4):
    with Pool(processes=n_processes) as pool:
        worker = partial(clusters_asm_par, cl=cl, fq_path=fq_path, output_dir=output_dir, flye_path=flye_path)
        list(tqdm(pool.imap(worker, ov_nodes), total=len(ov_nodes)))



def clusters_asm(cluster, cl):
    flye_dir= f"{StRainyArgs().output_intermediate}/flye_outputs/asm/{cluster}/asm"
    cluster_fastq=f"{StRainyArgs().output_intermediate}/flye_outputs/asm/{cluster}/cluster_reads.fastq"
    try:
        os.mkdir(f"{StRainyArgs().output_intermediate}/flye_outputs/asm/{cluster}")
    except FileExistsError:
        pass
    out_file = open(cluster_fastq, "w+")
    extracted_reads=cl.loc[cl["Cluster"] == cluster]["ReadName"].to_numpy()  # store read names
    with gzip.open(StRainyArgs().fq, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if record.id in extracted_reads:
                SeqIO.write(record, out_file, "fastq")
    out_file.close()
    flye = "/Users/ekaterina.kazantseva/strainy2/strainy/submodules/Flye/bin/flye"
    flye_cmd = f"{flye} " f" --pacbio-hifi {cluster_fastq}"  f" -o {flye_dir}"
    subprocess.check_output(flye_cmd, shell=True, capture_output=False, stderr=open(os.devnull, "w"))




def clusters_vis_stats(G, cl, I):
    """Creates connection graph vis and statistics"""
    cl.loc[cl['Cluster'] == 'NA', 'Cluster'] = 0
    cmap = plt.get_cmap('viridis')
    clusters=sorted(set(cl['Cluster'].astype(int)))
    cmap = cmap(np.linspace(0, 1, len(clusters)))
    colors = {}
    try:
        clusters.remove('0')
    except (KeyError, ValueError):
        pass
    colors[0] = "#505050"
    i = 0

    for cluster in clusters:
        colors[cluster] = mt.colors.to_hex(cmap[i])
        i = i + 1

    for index in cl.index:
        cl.loc[index, 'Color'] = colors[int(cl.loc[index, 'Cluster'])]
        G.remove_edges_from(list(nx.selfloop_edges(G)))
    try:
        nx.draw(G, nodelist=G.nodes(), with_labels=True, width=0.03, node_size=10, font_size=10,node_color=cl['Color'])
    except AttributeError:  #incompatability with scipy < 1.8
        pass




def color_bam2(edges):
    for edge in edges:
        color(edge, cl_file="%s/clusters/clusters_before_splitting_%s_%s.csv" %
                            (StRainyArgs().output_intermediate, I, StRainyArgs().AF),
              file="%s/bam/merged/coloredBAM_unitig_%s_merged.bam" % (StRainyArgs().output_intermediate, edge))
        out_bam_dir = os.path.join(StRainyArgs().output_intermediate, "bam/merged")
        final_aln = os.path.join(StRainyArgs().output, "alignment_phased_merged.bam")


    files_to_be_merged = []
    for fname in subprocess.check_output(f'find {out_bam_dir} -name "*unitig*.bam"',
                                         shell = True, universal_newlines = True).split("\n"):
        if len(fname):
            files_to_be_merged.append(fname)

    # Number of file to be merged could be > 4092,
    # in which case samtools merge throws too many open files error
    for i, bam_file in enumerate(files_to_be_merged):
        # fetch the header and put it at the top of the file, for the first bam_file only
        if i == 0:
            subprocess.check_output(f'samtools view -H {bam_file} > '
                                    f'{out_bam_dir}/coloredSAM.sam',shell = True)

        # convert bam to sam, append to the file
        subprocess.check_output(f'samtools view {bam_file} >> {out_bam_dir}/coloredSAM.sam',
                                shell = True)

    # convert the file to bam and sort
    subprocess.check_output(f'samtools view -b {out_bam_dir}/coloredSAM.sam >> '
                            f'{out_bam_dir}/unsortedBAM.bam',shell = True)
    pysam.samtools.sort(f'{out_bam_dir}/unsortedBAM.bam', "-o", final_aln)
    pysam.samtools.index(final_aln)

    # remove unnecessary files
    os.remove(f'{out_bam_dir}/unsortedBAM.bam')
    os.remove(f'{out_bam_dir}/coloredSAM.sam')
    for file in files_to_be_merged:
        os.remove(file)


if __name__ == "__main__":
    all_main()

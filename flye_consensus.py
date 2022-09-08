import subprocess
import os
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from params import *

"""
extract_reads
Created by Tim Stuart
based on https://timoast.github.io/blog/2015-10-12-extractreads/
"""
def extract_reads(read_names, bam_file, output_file, edge=""):
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
    header = bamfile.header.copy()
    cluster_start = -1
    cluster_end = -1
    read_list = []  # stores to reads to be written after the cluster start/end is calculated
    for name in read_names:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                if x.reference_name == edge or edge == "":
                    if x.reference_start < cluster_start or cluster_start == -1:
                        cluster_start = x.reference_start
                    if x.reference_end > cluster_end or cluster_end == -1:
                        cluster_end = x.reference_end
                    edge = x.reference_name
                    read_list.append(x)

    out = pysam.Samfile(output_file, "wb", header=header)
    for x in read_list:
        temp_dict = x.to_dict()
        temp_dict["ref_pos"] = str(int(temp_dict["ref_pos"]) - cluster_start)
        y = x.from_dict(temp_dict, x.header)  # create a new read from the modified dictionary
        out.write(y)

    out.close()
    return edge, cluster_start, cluster_end


# compute the consensus for a single cluster and return it as a string
def flye_consensus(edge, cl, cluster, data):
    # clusters: ID (int of the clusters)
    # cl: a .csv file under the clusters dir with fields (read name, cluster id)

    reads_from_curr_cluster = cl.loc[cl["Cluster"] == cluster]["ReadName"].to_numpy()  # store read names
    # TODO: use edge parameter
    edge, cluster_start, cluster_end = extract_reads(reads_from_curr_cluster, bam, f"cluster_{cluster}_reads.bam")
    g = gfapy.Gfa.from_file(gfa)  # read the gfa file
    # access the edge in the graph and cut according to the cluster start and end
    edge_seq_cut = (g.line(edge)).sequence[cluster_start:cluster_end]
    print(f"CLUSTER:{cluster}, CLUSTER_START:{cluster_start}, CLUSTER_END:{cluster_end}")
    # create a new fasta file with the sequence of the cut edge (will be an input to the Flye polisher)
    fname = f"{edge}-cluster{cluster}"
    record = SeqRecord(
        Seq(edge_seq_cut),
        id=f"{edge}",
        name=f"{edge} sequence cut for cluster {cluster}",
        description=""
    )
    SeqIO.write([record], f"output/{fname}.fa", "fasta")

    # sort the bam file
    sort_cmd = f"samtools sort cluster_{cluster}_reads.bam > cluster_{cluster}_reads_sorted.bam"
    subprocess.check_output(sort_cmd, shell=True, capture_output=False)
    # index the bam file
    index_cmd = f"samtools index cluster_{cluster}_reads_sorted.bam"
    subprocess.check_output(index_cmd, shell=True, capture_output=False)
    # run the Flye polisher
    polish_cmd = f"{flye} --polish-target output/{fname}.fa --pacbio-hifi cluster_{cluster}_reads_sorted.bam " \
                 f"-o output/flye_consensus_{edge}_{cluster}"
    subprocess.check_output(polish_cmd, shell=True, capture_output=False)

    # read back the output of the Flye polisher
    consensus = SeqIO.read(f"output/flye_consensus_{edge}_{cluster}/polished_1.fasta", "fasta")

    # delete the created input files for Flye
    os.remove(f"output/{fname}.fa")
    os.remove(f"cluster_{cluster}_reads.bam")
    os.remove(f"cluster_{cluster}_reads_sorted.bam")
    os.remove(f"cluster_{cluster}_reads_sorted.bam.bai")
    shutil.rmtree(f"output/flye_consensus_{edge}_{cluster}")

    return consensus.seq

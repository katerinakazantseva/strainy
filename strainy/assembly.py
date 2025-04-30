import gfapy
from strainy.params import StRainyArgs, init_global_args_storage
import strainy.graph_operations.gfa_ops as gfa_ops
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
from Bio import SeqIO
import subprocess
import gzip
import os


class StrainyAssembly:
    def __init__(self,overlap,cl):
        run_clusters_asm_parallel(
            overlap.nodes,
            cl,
            fq_path=StRainyArgs().fq,
            output_dir=StRainyArgs().output_intermediate,
            flye_path="/Users/ekaterina.kazantseva/strainy2/strainy/submodules/Flye/bin/flye",
            n_processes=StRainyArgs().threads
        )
        self.nodes = overlap.nodes
        self.edges = overlap.edges
        graph = gfapy.Gfa()
        for i in self.nodes:
            for record in SeqIO.parse(f"{StRainyArgs().output_intermediate}/flye_outputs/asm/{i}/asm/assembly.fasta",
                                  "fasta"):
                seq = record.seq
            gfa_ops.add_edge(graph, str(i), 1, seq)
        for i in self.edges:
            gfa_ops.add_link(graph, str(i[0]), "+", str(i[1]), "+", 1)
        self.graph = graph
    def write(self, path):
        gfapy.Gfa.to_file(self.graph,path)

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



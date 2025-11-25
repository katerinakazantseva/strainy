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
    #ADDD clipping
    fasta = f"{output_dir}/flye_outputs/asm/{cluster}/asm/assembly.fasta"
    bedgz = f"{output_dir}/flye_outputs/asm/{cluster}/asm/40-polishing/base_coverage.bed.gz"
    outfa = f"{output_dir}/flye_outputs/asm/{cluster}/asm/assembly_trimmed.fasta"
    threshold=3
    if not os.path.exists(fasta):
        print(f"[!] FASTA not found: {fasta}")
        return
    if not os.path.exists(bedgz):
        print(f"[!] BED.gz not found: {bedgz}")
        return
    p1 = subprocess.run(["gunzip", "-c", bedgz],
                            stdout=subprocess.PIPE, text=True, check=True)
    p2 = subprocess.run(["awk", "-v", f"t={threshold}", '$4 >= t {print $1"\t"$2"\t"$3}'],
                            input=p1.stdout, stdout=subprocess.PIPE, text=True, check=True)
    if not p2.stdout.strip():
        open(outfa, "w").close()
    p3 = subprocess.run(["bedtools", "merge", "-i", "-"],
                        input=p2.stdout, stdout=subprocess.PIPE, text=True, check=True)
    if not p3.stdout.strip():
        open(outfa, "w").close()
    subprocess.run(["bedtools", "getfasta", "-fi", fasta, "-bed", "-", "-fo", outfa],
                       input=p3.stdout, text=True, check=True)




def run_clusters_asm_parallel(ov_nodes, cl, fq_path, output_dir, flye_path, n_processes):
    with Pool(processes=n_processes) as pool:
        worker = partial(clusters_asm_par, cl=cl, fq_path=fq_path, output_dir=output_dir, flye_path=flye_path)
        list(tqdm(pool.imap(worker, ov_nodes), total=len(ov_nodes)))




class StrainyAssembly:
    def __init__(self,cl):
        nodes=set(cl["Cluster"].values)
        run_clusters_asm_parallel(
            nodes,
            cl,
            fq_path=StRainyArgs().fq,
            output_dir=StRainyArgs().output_intermediate,
            flye_path="/Users/katya/project_strainy/strainy2/strainy/submodules/Flye/bin/flye",
            n_processes=8
        )
        ofile = f"{StRainyArgs().output_intermediate}/contigs.fasta"
        ofile_trim = f"{StRainyArgs().output_intermediate}/contigs_trimmed.fasta"
        outfile = open(ofile, "w+")
        outfile_trimmed = open(ofile_trim, "w+")

        graph = gfapy.Gfa()
        for i in nodes:
            try:
                file = f"{StRainyArgs().output_intermediate}/flye_outputs/asm/{i}/asm/assembly.fasta"
                for record in SeqIO.parse(file, "fasta"):
                    new_id="cluster" + str(i)+str(record.id)
                    record.id = new_id
                    record.name = new_id
                    seq = record.seq
                    try:
                        gfa_ops.add_edge(graph, "cluster" + str(i)+str(record.description), 1, seq)
                    except:
                        print("ERROR")
                    SeqIO.write(record, outfile, "fasta")
            except FileNotFoundError:
                print("not found")
                continue
        outfile.close()
        for i in nodes:
            try:
                file = f"{StRainyArgs().output_intermediate}/flye_outputs/asm/{i}/asm/assembly_trimmed.fasta"
                for record in SeqIO.parse(file, "fasta"):
                    new_id="cluster" + str(i)+str(record.id)
                    record.id = new_id
                    record.name = new_id
                    seq = record.seq
                    try:
                        gfa_ops.add_edge(graph, "cluster" + str(i)+str(record.description), 1, seq)
                    except:
                        print("ERROR")
                    SeqIO.write(record, outfile_trimmed, "fasta")
            except FileNotFoundError:
                print("not found")
                continue
        outfile_trimmed.close()
        self.graph = graph


    def write(self, path):
        gfapy.Gfa.to_file(self.graph,path)

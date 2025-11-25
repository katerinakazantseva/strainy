import subprocess
import pysam
import os
import io
import re
from collections import Counter, namedtuple
from Bio import SeqIO
from strainy.params import *
import logging
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
import logging
import pandas as pd
import subprocess
import pysam
import os
import io
import re
from collections import Counter, namedtuple
from Bio import SeqIO
from strainy.params import *
import logging
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
import logging
import pandas as pd
logger = logging.getLogger()


def clusters(data):
    """
    Creates file for mapping ReadNames-Cluster N
    """
    cl = pd.DataFrame(columns=['ReadName', 'Cluster', 'Coordinates'])
    for key, value in data.items():
        coord = {k: [v["Start"], v["End"]] for k, v in value.items()}
        row = pd.DataFrame({'ReadName': [key], 'Cluster': ['NA'], 'Coordinates': [coord]})
        cl = pd.concat([cl, row])
    cl = cl.reset_index(drop=True)
    return cl

def read_snp2(vcf_file, bam, AF, cluster=None):
    """
       Extracts SNP positions from a VCF file or generates them from BAM data if no VCF file is provided.
       This function either reads SNP positions directly from a VCF file or generates SNP data from the BAM file
       using `bcftools` when no VCF file is provided. It filters SNPs based on allele frequency and read counts
       and writes filtered SNP information to a VCF file.
       Returns:
           list: A list of SNP positions that meet the filtering criteria.
       Notes:
           - If `vcf_file` is `None`, the function generates SNP data using `bcftools mpileup` and `bcftools query`
             commands, saving the results in a temporary file.
           - The SNPs are filtered based on read counts and allele frequency (`AF`), and the filtered SNP positions
             are added to the `snp_pos` list.
           - For each SNP, if both the forward and reverse allele frequencies are greater than or equal to 60% of
             the allele frequency (`AF`) threshold and have more than two supporting reads, the SNP is considered.
           - If a `vcf_file` is provided, SNP positions are read directly from this file using `bcftools`.
       """
    snp_pos = {}
    if vcf_file == None:
        if cluster == None:
            snpos = ('bcftools mpileup  {} --no-reference -I --no-version --annotate FORMAT/AD --annotate FORMAT/ADR --annotate FORMAT/ADF   2>/dev/null | bcftools query -f  "%CHROM %POS [ %AD %DP %ADR %ADF  %REF %ALT]\n"  >{}/vcf/vcf.txt').format(bam, StRainyArgs().output_intermediate)

            subprocess.check_output(snpos, shell=True, capture_output=False)
            filtered_file='{}/vcf/vcf_filtered.vcf'.format(StRainyArgs().output_intermediate)
            if not os.path.exists(filtered_file):
                vcf_file_f = open(filtered_file, "a+")
                vcf_file_f.write("##fileformat=VCFv4.2\n")
                vcf_file_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            with open("%s/vcf/vcf.txt" % (StRainyArgs().output_intermediate)) as f:
                lines = f.readlines()
                for line in lines:
                    try:
                        snp_freq = int(str(line.split()[2]).split(',')[2])
                        pos_cov = int(line.split()[3])
                        min_snp_freq = max(unseparated_cluster_min_reads, AF * pos_cov)
                        if snp_freq >= min_snp_freq:
                            var_freqF=0
                            var_freqR=0
                            FreqF=[int(i) for i in (line.split()[4].split(','))]
                            FreqF_s=sorted(FreqF,reverse=True)
                            FreqR=[int(i) for i in (line.split()[5].split(','))]
                            FreqR_s=sorted(FreqR,reverse=True)

                            if FreqF.index(FreqF_s[0]) in [1,2] and FreqF.index(FreqF_s[1]) in [1,2] and FreqR.index(FreqR_s[0]) in [1,2] and FreqR.index(FreqR_s[1]) in [1,2]:
                                try:
                                    dpF=sum([int(i) for i in list(line.split()[4].split(',')) if int(i)>1])
                                    altF=int(line.split()[4].split(',')[2])
                                    var_freqF=altF/dpF
                                    dpR = sum([int(i) for i in list(line.split()[5].split(',')) if int(i)>1])
                                    altR=int(line.split()[5].split(',')[2])
                                    var_freqR=altR/dpR
                                except ZeroDivisionError:
                                    continue
                            if var_freqF>=(AF)*0.6 and var_freqR >= (AF) * 0.6 and altF>2 and altR>2:
                                try:
                                    snp_pos[line.split()[0]].append(line.split()[1])
                                except (KeyError,AttributeError):
                                    snp_pos[line.split()[0]] = [line.split()[1]]
                                try:
                                    vcf_file_f.write(str(line.split()[0])+"\t"+str(line.split()[1])+"\t.\t"+str(line.split()[6])+"\t"+str(line.split()[7])+"\t.\tPASS\t.\n")
                                except: pass

                    except(IndexError):
                        pass
            try:
                vcf_file_f.close()
            except: pass
        else:
            raise Exception("Shouldn't happen")
    else: ##TODO test it
        bcftools_cmd = f"bcftools view -f PASS -H {vcf_file} {edge} --types snps"
        bcf_proc = subprocess.Popen(bcftools_cmd, shell=True, stdout=subprocess.PIPE)
        for line in io.TextIOWrapper(bcf_proc.stdout, encoding="utf-8"):
            snp_pos.append(line.split()[1])
    return snp_pos



def read_bam_chunk(bam, edge_chunk, snp_pos, min_mapping_quality, min_base_quality, min_al_len, max_aln_error):
    bamfile = pysam.AlignmentFile(bam, "rb")
    data = {}
    ref_lengths = dict(zip(bamfile.references, bamfile.lengths))
    #duplicates=[]
    #all_reads=[]
    CIGAR_SOFT = 4
    CIGAR_HARD = 5
    for edge in edge_chunk:
        try:
            for read in bamfile.fetch(edge):
                clipping = False
                aln_len = read.reference_end - read.reference_start
                aln_divergence = 0
                edge_len = ref_lengths[edge]

                if read.has_tag("de"):
                    aln_divergence = read.get_tag("de")
                for (op, size) in read.cigartuples:
                    if op in [CIGAR_SOFT, CIGAR_HARD] and size > max_clipping:
                        clipping = True

                if read.mapping_quality < min_mapping_quality or aln_divergence > max_aln_error:
                    continue

                # only allow single read alignment per unitig
                # if read.query_name in all_reads:
                # duplicates.append(read.query_name)

                # all_reads.append(read.query_name)
                try:
                    if read.query_name in data[edge] and read.is_supplementary == True:
                        continue
                except KeyError:
                    pass

                ALN_GAP = 100
                if (not clipping and aln_len > min_al_len) or \
                        (min(read.reference_start, edge_len - read.reference_end) < start_end_gap):

                    # data[read.query_name]={}
                    try:
                        data[read.query_name][edge] = {}
                    except KeyError:
                        data[read.query_name] = {}
                        data[read.query_name][edge] = {}
                    data[read.query_name][edge]["Start"] = read.reference_start
                    data[read.query_name][edge]["End"] = read.reference_end
                    data[read.query_name][edge]["Rclip"] = []
                    data[read.query_name][edge]["Lclip"] = []

                    if read.has_tag("SA"):
                        strand = "+" if not read.is_reverse else "-"
                        suppl_aln = [_parse_cigar(read.query_name, edge, read.reference_start, strand,
                                                  read.cigarstring, read.mapping_quality, edge_len)]
                        suppl_aln += [_parse_sa(read.query_name, sa_str, ref_lengths)
                                      for sa_str in read.get_tag("SA").split(";") if sa_str]
                        suppl_aln.sort(key=lambda a: a.query_start)

                        good_connections = []
                        for a1, a2 in zip(suppl_aln[:-1], suppl_aln[1:]):
                            if abs(a1.query_end - a2.query_start) < ALN_GAP and \
                                    min(a1.reference_start, a1.reference_length - a1.reference_end) < ALN_GAP and \
                                    min(a2.reference_start, a2.reference_length - a2.reference_end) < ALN_GAP:
                                good_connections.append((a1, a2))

                        for (a1, a2) in good_connections:
                            if a1.reference_name == a2.reference_name:
                                continue
                            if a1.reference_name == read.reference_name:
                                if a1.strand == "+":
                                    data[read.query_name][edge]["Rclip"].append((a2.reference_name, a2.strand))
                                else:
                                    data[read.query_name][edge]["Lclip"].append(
                                        (a2.reference_name, _neg_strand(a2.strand)))
                            if a2.reference_name == read.reference_name:
                                if a2.strand == "+":
                                    data[read.query_name][edge]["Lclip"].append((a1.reference_name, a1.strand))
                                else:
                                    data[read.query_name][edge]["Rclip"].append(
                                        (a1.reference_name, _neg_strand(a1.strand)))
            try:
                for pos in snp_pos[edge]:
                    for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools',
                                                       min_base_quality=min_base_quality,
                                                       ignore_overlaps=False, min_mapping_quality=min_mapping_quality,
                                                       ignore_orphans=False, truncate=True):
                        for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                try:
                                    if int(pos) >= data[pileupread.alignment.query_name][edge]["Start"] and int(pos) <= \
                                            data[pileupread.alignment.query_name][edge]["End"]:
                                        data[pileupread.alignment.query_name][edge][pos] = \
                                        pileupread.alignment.query_sequence[pileupread.query_position]

                                except KeyError:
                                    continue
            except KeyError:
                continue
        except ValueError:
            continue

    bamfile.close()
    return data

def read_bam2_parallel(
    bam, edges, snp_pos,
    min_mapping_quality,
    min_base_quality,
    min_al_len,
    max_aln_error,
    n_processes=None
):
    if n_processes is None:
        n_processes = min(8, cpu_count())

    logging.info(f"Using {n_processes} parallel processes to read BAM")

    chunk_size = len(edges) // n_processes + 1
    edge_chunks = [edges[i:i + chunk_size] for i in range(0, len(edges), chunk_size)]

    func = partial(
        read_bam_chunk,
        bam,
        snp_pos=snp_pos,
        min_mapping_quality=min_mapping_quality,
        min_base_quality=min_base_quality,
        min_al_len=min_al_len,
        max_aln_error=max_aln_error
    )

    merged_data = {}

    with Pool(processes=n_processes) as pool:
        for partial_data in tqdm(pool.imap(func, edge_chunks), total=len(edge_chunks), desc="Processing BAM chunks"):
            for k, v in partial_data.items():
                if k not in merged_data:
                    merged_data[k] = v
                else:
                    merged_data[k].update(v)

    logging.info("Finished reading BAM and extracting data.")
    return merged_data



def read_bam2(bam, edges, snp_pos, min_mapping_quality,min_base_quality, min_al_len, max_aln_error):
    """
     Extracts read alignment information from a BAM file for a specific edge, focusing on high-quality reads
     and their supplementary alignments.
     This function processes reads from a BAM file to gather information about their alignments on a given `edge`.
     It filters reads based on mapping quality, alignment length, and divergence, and collects details about
     clipping and supplementary alignments. Additionally, it captures base information at specified SNP positions.

     Returns:
         dict: A dictionary containing read alignment data for the specified `edge`. Each key is a read name,
               and its value is another dictionary with keys:
               - "Start": Start position of the read alignment.
               - "End": End position of the read alignment.
               - "Rclip": List of right-side clipping information for supplementary alignments.
               - "Lclip": List of left-side clipping information for supplementary alignments.
               - SNP positions as keys with corresponding base values as the read sequence at that position.
     """
    bamfile = pysam.AlignmentFile(bam, "rb")
    duplicates=[]
    all_reads=[]
    data = {}
    ref_lengths = dict(zip(bamfile.references, bamfile.lengths))

    CIGAR_SOFT = 4
    CIGAR_HARD = 5
    for edge in edges:
        try:
            for read in bamfile.fetch(edge):
                clipping = False
                aln_len = read.reference_end - read.reference_start
                aln_divergence = 0
                edge_len = ref_lengths[edge]

                if read.has_tag("de"):
                    aln_divergence = read.get_tag("de")
                for (op, size) in read.cigartuples:
                    if op in [CIGAR_SOFT, CIGAR_HARD] and size > max_clipping:
                        clipping = True

                if read.mapping_quality < min_mapping_quality or aln_divergence > max_aln_error:
                    continue


            #only allow single read alignment per unitig
            #if read.query_name in all_reads:
                #duplicates.append(read.query_name)

            #all_reads.append(read.query_name)
                try:
                    if read.query_name in data[edge] and read.is_supplementary==True:
                        continue
                except KeyError:
                    pass

                ALN_GAP = 100
                if (not clipping and aln_len > min_al_len) or \
                    (min(read.reference_start, edge_len - read.reference_end) < start_end_gap):

                    #data[read.query_name]={}
                    try:
                        data[read.query_name][edge] = {}
                    except KeyError:
                        data[read.query_name] = {}
                        data[read.query_name][edge] = {}
                    data[read.query_name][edge]["Start"] = read.reference_start
                    data[read.query_name][edge]["End"] = read.reference_end
                    data[read.query_name][edge]["Rclip"] = []
                    data[read.query_name][edge]["Lclip"] = []

                    if read.has_tag("SA"):
                        strand = "+" if not read.is_reverse else "-"
                        suppl_aln = [_parse_cigar(read.query_name, edge, read.reference_start, strand,
                                          read.cigarstring, read.mapping_quality, edge_len)]
                        suppl_aln += [_parse_sa(read.query_name, sa_str, ref_lengths)
                              for sa_str in read.get_tag("SA").split(";") if sa_str]
                        suppl_aln.sort(key=lambda a: a.query_start)

                        good_connections = []
                        for a1, a2 in zip(suppl_aln[:-1], suppl_aln[1:]):
                            if abs(a1.query_end - a2.query_start) < ALN_GAP and \
                                    min(a1.reference_start, a1.reference_length - a1.reference_end) < ALN_GAP and \
                                    min(a2.reference_start, a2.reference_length - a2.reference_end) < ALN_GAP:
                                good_connections.append((a1, a2))

                        for (a1, a2) in good_connections:
                            if a1.reference_name == a2.reference_name:
                                continue
                            if a1.reference_name == read.reference_name:
                                if a1.strand == "+":
                                    data[read.query_name][edge]["Rclip"].append((a2.reference_name, a2.strand))
                                else:
                                    data[read.query_name][edge]["Lclip"].append((a2.reference_name, _neg_strand(a2.strand)))
                            if a2.reference_name == read.reference_name:
                                if a2.strand == "+":
                                    data[read.query_name][edge]["Lclip"].append((a1.reference_name, a1.strand))
                                else:
                                    data[read.query_name][edge]["Rclip"].append((a1.reference_name, _neg_strand(a1.strand)))
            try:
                for pos in snp_pos[edge]:
                    for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools', min_base_quality=min_base_quality,
                                            ignore_overlaps=False, min_mapping_quality=min_mapping_quality,
                                            ignore_orphans=False, truncate=True):
                        for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                try:
                                    if int(pos) >= data[pileupread.alignment.query_name][edge]["Start"] and int(pos) <= data[pileupread.alignment.query_name][edge]["End"]:
                                        data[pileupread.alignment.query_name][edge][pos] = pileupread.alignment.query_sequence[pileupread.query_position]

                                except KeyError:
                                    continue
            except KeyError:
                continue
        except ValueError:
            continue
    bamfile.close()
    return data


def read_fasta_seq(filename, seq_name):
    reference_seq = None
    for seq in SeqIO.parse(filename, "fasta"):
        if seq.id == seq_name:
            reference_seq = str(seq.seq)
            break
    if reference_seq is None:
        raise Exception("Reference sequence not found")

    return reference_seq


def cluster_bounds(cl, cluster,edge):
    #TODO take into accout coverage (starting from two reads)
    coordinates = cl.loc[cl["Cluster"] == cluster]["Coordinates"].to_numpy()
    start_pos_of_reads = []
    end_pos_of_reads = []
    for i, c in enumerate(coordinates):
        try:
            start_pos_of_reads.append(c[edge][0])
            end_pos_of_reads.append(c[edge][1])
        except KeyError:
            continue
    try:
        bounds = [min(start_pos_of_reads),max(end_pos_of_reads)]
    except ValueError:
        bounds=[-1,-1]
    """
    bounds={}
    bounds[cluster] = {}
    val = {}
    for edge in edges:
        starts=[]
        ends=[]
        for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
            try:
                start = int(data[read][edge]["Start"])
                stop = int(data[read][edge]["End"])
                starts.append(start)
                ends.append(stop)
                print()
                #clCov = clCov + (stop - start)
            except(KeyError):
                continue
        try:
            clStart = sorted(starts)[1]
            clStop = sorted(ends)[len(ends) - 2]
            val["End"] = clStop
            val["Start"] = clStart
            bounds[cluster][edge] = val
        except(IndexError):
            pass

    """
    return bounds


def cluster_edges(cl, cluster):
    coordinates = cl.loc[cl["Cluster"] == cluster]["Coordinates"].to_numpy()
    list=[]
    for i,c in enumerate(coordinates):
        for k,v in c.items():
            list.append(k)
    return set(list)

cigar_parser = re.compile("[0-9]+[MIDNSHP=X]")
ReadSegment = namedtuple("ReadSegment", ["query_start", "query_end", "reference_start", "reference_end", "query_name", "reference_name",
                                         "strand", "reference_length", "query_length", "mapq"])




def _parse_cigar(read_id, ref_id, ref_start, strand, cigar, mapq, ref_length):
    """
    Parses cigar and generate ReadSegment structure with alignment coordinates
    """
    first_clip = True
    read_start = 0
    read_aligned = 0
    read_length = 0
    ref_aligned = 0
    ref_start = int(ref_start)
    mapq = int(mapq)

    for token in cigar_parser.findall(cigar):
        op = token[-1]
        op_len = int(token[:-1])

        if op == "H" or op == "S":
            if first_clip:
                read_start = op_len
            read_length += op_len
        first_clip = False

        if op == "M" or op == "=" or op == "X":
            read_aligned += op_len
            ref_aligned += op_len
            read_length += op_len
        if op == "D":
            ref_aligned += op_len
        if op == "I":
            read_aligned += op_len
            read_length += op_len

    ref_end = ref_start + ref_aligned
    read_end = read_start + read_aligned

    if strand == "-":
        read_start, read_end = read_length - read_end, read_length - read_start

    return ReadSegment(read_start, read_end, ref_start, ref_end, read_id,
                       ref_id, strand, ref_length, read_length, mapq)


def _parse_sa(read_id, sa_str, ref_lengths):
    ref_id, ref_start, strand, cigar, mapq, _nm = sa_str.split(",")
    return _parse_cigar(read_id, ref_id, ref_start, strand, cigar, mapq, ref_lengths[ref_id])




def _neg_strand(strand):
    return "-" if strand == "+" else "+"



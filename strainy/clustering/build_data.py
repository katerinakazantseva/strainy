import subprocess
import pysam
import os
import io
import re
from collections import Counter, namedtuple
from Bio import SeqIO
from strainy.params import *
import logging
logger = logging.getLogger()




def read_snp(vcf_file, edge, bam, AF, cluster=None):
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
    snp_pos = []
    if vcf_file == None:
        if cluster == None:
            snpos = ('bcftools mpileup -r {} {} --no-reference -I --no-version --annotate FORMAT/AD --annotate FORMAT/ADR --annotate FORMAT/ADF   2>/dev/null | bcftools query -f  "%CHROM %POS [ %AD %DP %ADR %ADF  %REF %ALT]\n"  >{}/vcf/vcf_{}.txt').format(edge, bam, StRainyArgs().output_intermediate, edge)

            subprocess.check_output(snpos, shell=True, capture_output=False)
            filtered_file='{}/vcf/vcf_{}_filtered.vcf'.format(StRainyArgs().output_intermediate, edge)
            if not os.path.exists(filtered_file):
                vcf_file_f = open(filtered_file, "a+")
                vcf_file_f.write("##fileformat=VCFv4.2\n")
                vcf_file_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            with open("%s/vcf/vcf_%s.txt" % (StRainyArgs().output_intermediate, edge)) as f:
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
                                snp_pos.append(line.split()[1])
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
    else:
        bcftools_cmd = f"bcftools view -f PASS -H {vcf_file} {edge} --types snps"
        bcf_proc = subprocess.Popen(bcftools_cmd, shell=True, stdout=subprocess.PIPE)
        for line in io.TextIOWrapper(bcf_proc.stdout, encoding="utf-8"):
            snp_pos.append(line.split()[1])
    return snp_pos


ReadSegment = namedtuple("ReadSegment", ["query_start", "query_end", "reference_start", "reference_end", "query_name", "reference_name",
                                         "strand", "reference_length", "query_length", "mapq"])
cigar_parser = re.compile("[0-9]+[MIDNSHP=X]")




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




def read_bam(bam, edge, snp_pos, min_mapping_quality,min_base_quality, min_al_len, max_aln_error):
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
        if read.query_name in data and read.is_supplementary==True:
            continue

        ALN_GAP = 100
        if (not clipping and aln_len > min_al_len) or \
                (min(read.reference_start, edge_len - read.reference_end) < start_end_gap):
            data[read.query_name] = {}
            data[read.query_name]["Start"] = read.reference_start
            data[read.query_name]["End"] = read.reference_end
            data[read.query_name]["Rclip"] = []
            data[read.query_name]["Lclip"] = []

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
                            data[read.query_name]["Rclip"].append((a2.reference_name, a2.strand))
                        else:
                            data[read.query_name]["Lclip"].append((a2.reference_name, _neg_strand(a2.strand)))
                    if a2.reference_name == read.reference_name:
                        if a2.strand == "+":
                            data[read.query_name]["Lclip"].append((a1.reference_name, a1.strand))
                        else:
                            data[read.query_name]["Rclip"].append((a1.reference_name, _neg_strand(a1.strand)))

    for pos in snp_pos:
        for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools', min_base_quality=min_base_quality,
                                           ignore_overlaps=False, min_mapping_quality=min_mapping_quality,
                                           ignore_orphans=False, truncate=True):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    try:
                        if int(pos) >= data[pileupread.alignment.query_name]["Start"] and int(pos) <= data[pileupread.alignment.query_name]["End"]:
                            data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[pileupread.query_position]

                    except (KeyError):
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




def build_data_cons(cl, snp_pos, data, edge, reference_seq):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
    cons = {}
    for cluster in clusters:
        cons = cluster_consensuns(cl, cluster, snp_pos, data, cons, edge, reference_seq)
    return cons




def cluster_consensuns(cl, cluster, snp_pos, data, cons, edge, reference_seq):
    """
    Generates a consensus sequence for a given cluster based on SNP positions and read data.
    This function processes a cluster of reads, extracts SNP positions, and determines the most frequent bases
    (alleles) to form a consensus sequence. It also identifies potential variants and calculates statistics such as
    coverage, start, and stop positions of the cluster. The results are stored in the `cons` dictionary.
    Returns:
        dict: The updated `cons` dictionary containing consensus information for the specified cluster.
    """
    strange = 0
    strange2 = 0
    val = {}
    clust_snp = []
    mpileup_snps = []
    mis_count=0
    Rcl=StRainyArgs().Rcl
    AF=StRainyArgs().AF
    for pos in snp_pos:
        npos = []
        for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
            try:
                npos.append(data[read][pos])
            except(KeyError):
                continue

        min_snp_freq = max(unseparated_cluster_min_reads, AF * len(npos))
        alt_snp_freq = max(unseparated_cluster_min_reads, split_allele_freq * len(npos))
        try:
            if len(npos) >= unseparated_cluster_min_reads:
                #store most frequent symbol as consensus
                if int(Counter(npos).most_common()[0][1]) > 2:
                    val[pos] = Counter(npos).most_common()[0][0]

                #mimicking bcftools mpileup
                for elem, freq in Counter(npos).most_common():
                    if elem != reference_seq[int(pos) - 1] and freq >= min_snp_freq:
                        mpileup_snps.append(pos)
                        break

            #2nd most frequent, indicating a variant
            if int(Counter(npos).most_common()[1][1]) >= alt_snp_freq:
                #print(cluster, pos, Counter(npos).most_common())
                #strange = 1
                mis_count=mis_count+1
                clust_snp.append(pos)

        except IndexError:
            continue
    

    clust_snp2 = mpileup_snps
    val["clust_snp"] = clust_snp
    val["clust_snp2"] = clust_snp2

    clStart = 1000000000000  # change fo ln
    clStop = 0
    clCov = 0

    starts=[]
    ends=[]

    for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
        try:
            start=int(data[read]["Start"])
            stop=data[read]["End"]
            starts.append(start)
            ends.append(stop)
            clCov = clCov + (stop - start)
        except(KeyError):
            pass
    try:
        clStart = sorted(starts)[1]
        clStop = sorted(ends)[len(ends)-2]
    except(IndexError):
        pass
    try:
        if len(clust_snp2) > 0 and max([int(clust_snp2[i + 1]) - int(clust_snp2[i])
                                    for i in range(0, len(clust_snp2) - 1)]) > 1.5 * I:
            strange2 = 1

        if (int(clust_snp2[0]) - int(clStart)) > 1.5 * I or \
                int(clStop) - int(clust_snp2[len(clust_snp2) - 1]) > 1.5 * I:
            strange2 = 1

    except (ValueError, IndexError):
        pass


    try:
        if mis_count / (clStop - clStart) > Rcl:
            strange = 1
    except ZeroDivisionError:
        strange = 1

    val["Strange"] = int(strange == 1)
    val["Strange2"] = int(strange2 == 1)
    clCov = clCov / (clStop - clStart) if clStop > clStart else 0
    val["End"] = clStop
    val["Start"] = clStart
    val["Cov"] = clCov
    cons[cluster] = val

    return cons

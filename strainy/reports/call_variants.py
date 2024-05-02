#!/usr/bin/env python


import sys
import subprocess
import numpy as np
from datetime import date


def segment_match(match_str):
    match_list = []
    indel_list = []
    pos = 0

    while pos < len(match_str):
        if match_str[pos] in b"^":   #start of the read + mapq
            pos += 2
            continue

        if match_str[pos] in b"$":  #end of the read
            pos += 1
            continue

        if match_str[pos] in b".,ACGTNacgtn":  #match fwd / revserse
            match_list.append((b"M", match_str[pos : pos + 1]))
            pos += 1
            continue

        if match_str[pos] in b"<>*#":   #deletion/ref skip
            match_list.append((b"S", ""))
            pos += 1
            continue

        if match_str[pos] in b"+-": #INS/DEL sequnece
            num_pos = pos + 1
            while match_str[num_pos] in b"0123456789":
                num_pos += 1
            indel_size = int(match_str[pos + 1 : num_pos])
            indel_seq = match_str[num_pos : num_pos + indel_size]
            indel_list.append((len(match_list) - 1, match_str[pos : pos + 1], indel_seq))
            pos = num_pos + indel_size
            continue

        print(pos, match_str[pos], b"^")
        raise Exception("Parsing error")

    return match_list, indel_list


def _vcf_header(out_stream, ref_fasta, sample_name):
    out_stream.write("##fileformat=VCFv4.2\n")
    timestamp = date.today()
    out_stream.write(f"##fileDate={timestamp}\n")
    out_stream.write(f"##source=Strainy\n")
    out_stream.write(f"##reference=file://{ref_fasta}\n")
    out_stream.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
    out_stream.write("##INFO=<ID=ALT_HAP,Number=1,Type=String,Description=\"Haplotypes supporting ALT\">\n")
    out_stream.write("##INFO=<ID=REF_HAP,Number=1,Type=String,Description=\"Haplotypes supporting REF\">\n")
    out_stream.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    out_stream.write("##FORMAT=<ID=DV,Number=1,Type=String,Description=\"Varinat depth\">\n")
    out_stream.write("##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Total depth\">\n")
    out_stream.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")


SAMTOOLS = "samtools"
def generate_vcf(bam_aln, ref_fasta, out_stream):
    cmd = [SAMTOOLS, "mpileup", bam_aln, "-f", ref_fasta, "--output-MQ", "--output-extra", "QNAME", "--min-MQ", "10"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    _vcf_header(out_stream, ref_fasta, "sample")

    strainy_id = 0
    for line in proc.stdout:
        fields = line.split()
        if len(fields) < 8:
            continue

        match_str = fields[4]
        if not(match_str) or match_str == b"." * len(match_str):
            continue

        ref_id, ref_pos, ref_base, num_reads, mapqs, read_names = \
                fields[0], int(fields[1]), fields[2], int(fields[3]), fields[6], fields[7]
        match_list, indel_list = segment_match(match_str)
        if len(match_list) != num_reads:
            raise Exception("Length mismatch")

        read_list = read_names.decode("utf-8").split(",")
        qual_converted = [int(c) - 33 for c in mapqs]
        qual = int(np.median(qual_converted))

        #INDEL
        if indel_list:
            ref_str = None
            alt_str = None

            read_support = []
            for indel in indel_list:
                read_support.append(read_list[indel[0]])
            vcf_dv = len(read_support)

            if indel_list[0][1] == b"+":
                vcf_id = f"Strainy_INS_{strainy_id}"
                ref_str = ref_base.decode("utf-8")
                alt_str = ref_str + indel_list[0][2].decode("utf-8")
            else:
                vcf_id = f"Strainy_DEL_{strainy_id}"
                alt_str = ref_base.decode("utf-8")
                ref_str = alt_str + indel_list[0][2].decode("utf-8")

            info = "ALT_HAP=" + ",".join(read_support)
            ref_hap = list(set(read_list) - set(read_support))
            if ref_hap:
                info += ";REF_HAP=" + ",".join(ref_hap)

            vcf_dp = len(read_list)
            vcf_fields = [ref_id.decode("utf-8"), ref_pos, vcf_id, ref_str, alt_str,
                          qual, "PASS", info, "GT:DV:DP", f"0/1:{vcf_dv}:{vcf_dp}"]
            out_stream.write("\t".join(map(str, vcf_fields)) + "\n")

            strainy_id += 1

        #SNP
        else:
            read_support = []
            alt_snps = set()
            for i, read in enumerate(match_list):
                if read[0] in b"M" and read[1] in b"ACGTacgt":
                    alt_snps.add(read[1].decode("utf-8"))
                    read_support.append(read_list[i])

            if len(alt_snps):
            #if len(alt_snps) and len(read_support) > 1:
                vcf_id = f"Strainy_SNP_{strainy_id}"
                strainy_id += 1

                vcf_dp = len(read_list)
                vcf_dv = len(read_support)
                info = "ALT_HAP=" + ",".join(read_support)
                ref_hap = list(set(read_list) - set(read_support))
                if ref_hap:
                    info += ";REF_HAP=" + ",".join(ref_hap)

                vcf_fields = [ref_id.decode("utf-8"), ref_pos, vcf_id, ref_base.decode("utf-8"), ",".join(list(alt_snps)),
                              qual, "PASS", info, "GT:DV:DP", f"0/1:{vcf_dv}:{vcf_dp}"]
                out_stream.write("\t".join(map(str, vcf_fields)) + "\n")


MINIMAP = "minimap2"
def run_minimap2(reference, query, threads, aln_out):
    cmd = f"{MINIMAP} -ax asm10 {reference} {query} -t {threads} | {SAMTOOLS} sort -@4 > {aln_out}"
    subprocess.check_call(cmd, shell=True)
    subprocess.check_call(f"{SAMTOOLS} index -@4 {aln_out}", shell=True)


def produce_strainy_vcf(ref_fasta, strain_asm, threads, minimap_aln, vcf_stream):
    run_minimap2(ref_fasta, strain_asm, threads, minimap_aln)
    generate_vcf(minimap_aln, ref_fasta, vcf_stream)

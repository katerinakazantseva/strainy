#!/usr/bin/env python3

import sys
import csv
import subprocess
from collections import defaultdict, namedtuple
import math

import numpy as np


RefInfo = namedtuple("RefInfo", ["ref_utg", "length", "snp_rate", "is_processed", "is_phased"])
PhasedInfo = namedtuple("PhasedInfo", ["ref_utg", "pahsed_utg", "phased_len", "cov", "snp_rate", "ref_start", "ref_end"])

WND = 100
MIN_DIV = 0.00001


def _calc_n50(scaffolds_lengths, assembly_len):
    n50 = 0
    sum_len = 0
    l50 = 0
    for l in sorted(scaffolds_lengths, reverse=True):
        sum_len += l
        l50 += 1
        if sum_len > assembly_len / 2:
            n50 = l
            break
    return l50, n50


def get_lengths_and_coverage(ref_utgs_table):
    ref_utgs = {}

    for line in open(ref_utgs_table, "r"):
        if line.startswith("Reference"):
            continue
        fields = line.strip().split("\t")
        name, length, coverage, snp_rate, processed, phased = \
            fields[0], int(fields[1]), int(fields[2]), float(fields[3]), fields[4] == "True", fields[5] == "True"
        ref_utgs[name] = RefInfo(name, length, snp_rate, processed, phased)

    return ref_utgs


def stats_by_ref(strainy_stats, ref_utgs_info, ref_whitelist, out_stream):

    all_cov_frequencies = defaultdict(int)
    all_snp_rates = defaultdict(int)
    all_depths = defaultdict(int)

    strains_by_ref = defaultdict(list)
    #for line in csv.reader(open(strainy_stats, "r")):
    for line in open(strainy_stats, "r"):
        if line.startswith("Strain"):
            continue
        line = line.strip().split("\t")

        ref_utg, strain_utg, strain_len, cov, snp_rate, strain_start, strain_end = \
                line[1], line[0], int(line[2]), int(line[3]), float(line[6]), int(line[7]), int(line[8])
        #strain_len = strain_end - strain_start

        #print(ref_utg, strain_utg, strain_len, snp_rate)
        strains_by_ref[ref_utg].append(PhasedInfo(ref_utg, strain_utg, strain_len, cov, snp_rate, strain_start, strain_end))

    total_strain_utgs = []
    #print("#ref_utg\tref_len\tref_cov\tnum_stu\tmean_cov\tmed_snp\tmed_strain_depth")
    for ref_utg in sorted(ref_utgs_info, key=lambda x: ref_utgs_info[x].length, reverse=True):
        ref_len = ref_utgs_info[ref_utg].length
        strain_utgs = strains_by_ref[ref_utg]
        if len(strain_utgs) == 0:
            continue

        if ref_whitelist and ref_utg not in ref_whitelist:
            continue

        strain_lengths = [s.phased_len for s in strain_utgs]
        strain_snp_rates = [s.snp_rate for s in strain_utgs]
        strain_covs = [s.cov for s in strain_utgs]
        total_coverage = sum(strain_covs)
        total_strain_utgs.extend(strain_lengths)

        coverage_bins = [0 for _ in range(0, ref_len // WND)]
        for st in strain_utgs:
            for i in range(st.ref_start // WND, st.ref_end // WND):
                coverage_bins[i] += 1
        #print(ref_utg)
        #print(coverage_bins)

        for cov, ctg_len in zip(strain_covs, strain_lengths):
            all_depths[cov] += ctg_len // WND

        for x in coverage_bins:
            all_cov_frequencies[x] += 1

        for rate, ctg_len in zip(strain_snp_rates, strain_lengths):
            if rate < MIN_DIV:
                continue
            log_x = -10 * math.log(rate, 10)
            all_snp_rates[int(log_x)] += ctg_len // WND

        mean_strain_cov = int(sum(strain_lengths) / ref_len)
        median_strain_cov = int(np.median(coverage_bins))
        vals = [ref_utg, ref_len, total_coverage, len(strain_utgs), median_strain_cov,
               np.median(strain_snp_rates), int(np.median(strain_covs))]
        #for st in strain_utgs:
        #   print("\t", st[0], st[3], st[4])
        #print(coverage_bins)
        #print(strain_snp_rates)

        #print("\t".join(map(str, vals)))

    ######

    ctgs_input = []
    ctgs_processed = []
    ctgs_phased = []

    for utg, info in ref_utgs_info.items():
        if ref_whitelist and utg not in ref_whitelist:
            continue

        ctgs_input.append(info.length)
        if info.is_processed:
            ctgs_processed.append(info.length)
        if info.is_phased:
            ctgs_phased.append(info.length)

    def _n50(x):
        return _calc_n50(x, sum(x))[1]

    out_stream.write("Reference utgs input:\tlen: {0}\tnum: {1}\tN50:{2}\n"
                        .format(sum(ctgs_input), len(ctgs_input), _n50(ctgs_input)))
    out_stream.write("Reference utgs select:\tlen: {0}\tnum: {1}\tN50:{2}\n"
                        .format(sum(ctgs_processed), len(ctgs_processed), _n50(ctgs_processed)))
    out_stream.write("Reference utgs phased:\tlen: {0}\tnum: {1}\tN50:{2}\n"
                        .format(sum(ctgs_phased), len(ctgs_phased), _n50(ctgs_phased)))
    out_stream.write("Strain utgs asmembled:\tlen: {0}\tnum: {1}\tN50:{2}\n"
                        .format(sum(total_strain_utgs), len(total_strain_utgs), _n50(total_strain_utgs)))

    out_stream.write("\nMultiplicity\n")
    out_stream.write("Mul\tRefSeqLength\n")
    max_hist = max(all_cov_frequencies.values())
    for i in sorted(all_cov_frequencies):
        ref_len = all_cov_frequencies[i] * WND
        rate = all_cov_frequencies[i] / max_hist
        hist = int(rate * 20) * "*"
        #strain_len = all_cov_frequencies[i] * WND * i
        if i > 0 and i <= 10:
            out_stream.write(f"{i}\t{ref_len:10}\t{hist}\n")

    out_stream.write("\nRead depth\n")
    max_hist = max(all_depths.values())
    for i in sorted(all_depths):
        depth = all_depths[i] * WND
        rate = all_depths[i] / max_hist
        hist = int(rate * 20) * "*"
        if rate > 0.05:
            out_stream.write(f"{i}\t{depth:10}\t{hist}\n")

    out_stream.write("\nSNP divergence\n")
    out_stream.write("Rate\tStrainSeq\n")
    max_hist = max(all_snp_rates.values())
    for i in sorted(all_snp_rates):
        div_perc = math.pow(10, -i / 10)
        if div_perc > 1:
            continue
        strain_seq = all_snp_rates[i] * WND
        rate = all_snp_rates[i] / max_hist
        hist = int(rate * 20) * "*"
        if rate > 0.01:
            out_stream.write(f"{div_perc:4.5f}\t{strain_seq}\t{hist}\n")

    #####


def strain_stats_report(ref_utg_stats, phased_utg_stats, out_stream):
    ref_utg_stats = get_lengths_and_coverage(ref_utg_stats)
    stats_by_ref(phased_utg_stats, ref_utg_stats, None, out_stream)

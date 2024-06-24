import re
import gfapy
from collections import Counter, deque, defaultdict
import pandas as pd
import pysam
import csv
from strainy.params import *
from strainy.clustering import build_data




def write_phased_unitig_csv():
    columns = [
        'Strain_unitig',
        'Reference_unitig',
        'Length',
        'Coverage',
        'Abundance_Ratio',
        '#SNP',
        'SNP_density',
        'Start_positioin',
        'End_position'
    ]

    with open(StRainyArgs().phased_unitig_info_table_path, 'w') as f:
        write = csv.writer(f, delimiter='\t')

        write.writerow(columns)
        write.writerows(list(StRainyArgs().phased_unitig_info_table.values()))


def write_reference_unitig_csv():
    columns =[
        'Reference_unitig',
        'Length',
        'Coverage',
        'SNP_density',
        'Is_processed',
        'Is_phased'
    ]

    with open(StRainyArgs().reference_unitig_info_table_path, 'w') as f:
        write = csv.writer(f, delimiter='\t')

        write.writerow(columns)
        write.writerows(list(StRainyArgs().reference_unitig_info_table.values()))


def format_rounding(number):
    n = abs(number)
    if n == 0:
        return '0.000'
    if n < 1:
        # Find the first non-zero digit.
        # We want 3 digits, starting at that location.
        s = f'{n:.99f}'
        index = re.search('[1-9]', s).start()
        return s[:index + 3]

    # We want 2 digits after decimal point.
    return str(round(n, 2))


def store_phased_unitig_info(strain_unitig, reference_unitig, n_SNPs, start, end):
    reference_coverage = round(float(pysam.samtools.coverage("-r",
                                                             reference_unitig,
                                                             StRainyArgs().bam,
                                                             "--no-header").
                                     split()[6]))
    # # Log the information to std output
    # logger.info(f'== == Inserted Strain unitig: {strain_unitig.name} == == ')
    # logger.info(f'\t\t Reference unitig: {reference_unitig}')
    # logger.info(f'\t\t Length: {strain_unitig.length} bp')
    # logger.info(f'\t\t Coverage: {strain_unitig.dp}')
    # logger.info(f'\t\t Abundance Ratio: {round(100 * strain_unitig.dp // reference_coverage)}%')
    # logger.info(f'\t\t #SNP: {n_SNPs}')
    # logger.info(f'\t\t SNP density: {format_rounding(n_SNPs / strain_unitig.length)}')
    # logger.info(f'\t\t Start position: {start}')
    # logger.info(f'\t\t End position: {end}\n\n')

    try:
        abundance_ratio = round(100 * strain_unitig.dp // reference_coverage)
    except ZeroDivisionError:
        abundance_ratio = 0

    StRainyArgs().phased_unitig_info_table[strain_unitig.name] = [
        strain_unitig.name,
        reference_unitig,
        strain_unitig.length,
        strain_unitig.dp,
        abundance_ratio,
        n_SNPs,
        format_rounding(n_SNPs / strain_unitig.length),
        start,
        end
    ]


def store_reference_unitig_info(ref_coverage):
    graph = gfapy.Gfa.from_file(StRainyArgs().gfa)
    phased_unitig_df = pd.read_csv(StRainyArgs().phased_unitig_info_table_path, sep='\t')
    counter = Counter(list(phased_unitig_df['Reference_unitig']))
    for reference_unitig in graph.segments:
        # Number of phased unitigs created from this reference unitig
        n_phased_unitigs = counter[reference_unitig.name]
        # Number of SNPs
        n_SNPs = len(
            build_data.read_snp(
                StRainyArgs().snp,
                reference_unitig.name,
                StRainyArgs().bam,
                StRainyArgs().AF
            )
        )
        StRainyArgs().reference_unitig_info_table[reference_unitig.name] = [
            reference_unitig.name,
            reference_unitig.length,
            ref_coverage[reference_unitig.name],
            format_rounding(n_SNPs / reference_unitig.length),
            reference_unitig.name in StRainyArgs().edges_to_phase,
            n_phased_unitigs > 1
        ]
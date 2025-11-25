import csv
import os
import logging
import gfapy
import matplotlib.pyplot as plt
import matplotlib as mt
import pysam
import pandas as pd
import numpy as np
from strainy.params import *
import subprocess
logging.getLogger("matplotlib.font_manager").disabled = True


def write_bam(edge, cl, infile,outfile):
    """Creates new bam file based on ifnfile and add YC tag to the alignment based on csv file"""
    iterbam = infile.fetch(edge,until_eof=True)
    cmap = plt.get_cmap("viridis")
    cl.loc[cl["Cluster"] == "NA", "Cluster"] = 10000
    clusters = set(cl["Cluster"])
    cmap = cmap(np.linspace(0, 1, len(clusters)))
    colors={}
    i=0
    colors[10000] = "#505050"
    try:
        clusters.remove("10000")
    except KeyError:
        pass
    for cluster in clusters:
        colors[cluster] = mt.colors.to_hex(cmap[i])
        i = i+1
    cl_dict = dict(zip(cl.ReadName, cl.Cluster))
    for read in iterbam:
        try:
            cl_n = cl_dict[str(read).split()[0]]
            tag = colors[cl_n]
            read.set_tag("YC", tag, replace=False)
            outfile.write(read)
        except KeyError:
            continue
    outfile.close()


def color(edge,cl_file=None,file=None):
    #Creates colored edge bam based on strainy csv file with clusters IDs by default
    try:
        infile = pysam.AlignmentFile(StRainyArgs().bam, "rb")
        if file is None:
            outfile = pysam.AlignmentFile(
                f"{StRainyArgs().output_intermediate}/bam/coloredBAM_unitig_{edge}.bam",
                "wb", template=infile)
        else:
            outfile = pysam.AlignmentFile(file,"wb", template=infile)
        if cl_file is None:
            cl = pd.read_csv(
                f"{StRainyArgs().output_intermediate}/clusters/clusters_{edge}_{I}_{StRainyArgs().AF}.csv",
                keep_default_na=False)
        else:
            cl = pd.read_csv(cl_file,keep_default_na=False)
        write_bam(edge,cl,infile,outfile)
    except FileNotFoundError:
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

"""
def clusters_vis_stats(G, cl, I):
    #Creates connection graph vis and statistics
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


"""


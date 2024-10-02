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
logging.getLogger("matplotlib.font_manager").disabled = True


def write_bam(edge, cl, infile,outfile):
    """Creates new bam file based on ifnfile and add YC tag to the alignment based on csv file"""
    iterbam = infile.fetch(edge,until_eof=True)
    cmap = plt.get_cmap("viridis")
    cl.loc[cl["Cluster"] == "NA", "Cluster"] = 0
    clusters = set(cl["Cluster"])
    cmap = cmap(np.linspace(0, 1, len(clusters)))
    colors={}
    i=0
    colors[0] = "#505050"
    try:
        clusters.remove("0")
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
    """Creates colored edge bam based on strainy csv file with clusters IDs by default"""
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

import csv
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mt
from params import *
import gfapy
import os




def write_bam(edge, I, AF):
    infile = pysam.AlignmentFile(bam, "rb")
    outfile = pysam.AlignmentFile("output/bam/coloredBAM_%s.bam" % edge, "wb", template=infile)
    cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF),keep_default_na=False)
    iter = infile.fetch(edge,until_eof=True)
    cmap = plt.get_cmap('viridis')
    cl.loc[cl['Cluster'] == 'NA', 'Cluster'] = 0
    clusters=sorted(set(cl['Cluster'].astype(int)))
    cmap = cmap(np.linspace(0, 1, len(clusters)))
    colors={}
    i=0
    colors[0] = "#505050"
    try:
        clusters.remove('0')
    except: KeyError
    for cluster in clusters:
        colors[cluster]=mt.colors.to_hex(cmap[i])
        i=i+1
    cl_dict = dict(zip(cl.ReadName, cl.Cluster))
    for read in iter:
        try:
            clN=int(cl_dict[str(read).split()[0]])
            tag=colors[clN]
            read.set_tag('YC', tag, replace=False)
            outfile.write(read)
        except (KeyError):
            continue
    outfile.close()

'''def index(outfile,edge, I, AF):
    print("index")
    index = "output/bam/coloredBAM_%s_%s_%s.bai" % (edge, I, AF)
    pysam.samtools.index(outfile, index)'''




def color(i):
    edge=edges[i]

    try:
        write_bam(edge, I, AF)
        print("done")
    except (FileNotFoundError):
        pass



'''
for edge in edges:
    print(edge)
    try:
        write_bam(edge, I, AF)
        print("done")
    except (FileNotFoundError):
        continue
infile.close()
outfile.close()

pysam.samtools.index("output/bam/coloredBAM.bam","output/bam/coloredBAM.bai")'''
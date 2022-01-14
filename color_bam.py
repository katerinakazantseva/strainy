import csv
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mt



edge='edge_1'


infile = pysam.AlignmentFile("/Users/ekaterina.kazantseva/MT/test_data/test.bam", "rb")


def write_bam(infile, edge):
    outfile = pysam.AlignmentFile("output/bam/coloredBAM_%s.bam" % edge, "wb", template=infile)
    cl = pd.read_csv("output/clusters_%s.csv" % edge)
    iter = infile.fetch(edge,until_eof=True)
    cmap = plt.get_cmap('viridis')
    cl.dropna(inplace=True)
    cl.Cluster = cl.Cluster.astype(int)
    clusters=set(cl['Cluster'])
    cmap = cmap(np.linspace(0, 1, len(clusters)))
    colors={}
    i=0
    for cluster in clusters:
        colors[cluster]=mt.colors.to_hex(cmap[i])
        i=i+1
    cl_dict = dict(zip(cl.ReadName, cl.Cluster))

    for read in iter:
        try:
            clN=cl_dict[str(read).split()[0]]
            tag=colors[clN]
            read.set_tag('YC', tag, replace=False)
            outfile.write(read)
        except (KeyError):
            continue
    infile.close()
    outfile.close()

def index(outfile):
    index = "output/bam/coloredBAM_%s.bai" % edge
    pysam.samtools.index(outfile, index)




write_bam(infile, edge)
index("output/bam/coloredBAM_%s.bam" % edge)
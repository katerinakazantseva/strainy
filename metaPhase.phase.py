import multiprocessing
from params import *
from cluster import cluster
from color_bam import color
import pysam
import os
import subprocess




#dirs = ( "%s/vcf/" % output ,"output/adj_M/","output/clusters/","output/graphs/","output/bam/","output/bam/clusters")
dirs = ( "%s/vcf/" % output ,"%s/adj_M/" % output,"%s/clusters/" % output,"%s/graphs/" % output,"%s/bam/" % output,"%s/bam/clusters" % output)
for dir in dirs:
    try:
        os.stat(dir)
    except:
         os.makedirs(dir)


def phase(edges):
    pool = multiprocessing.Pool(3)
    pool.map(cluster, range(0, len(edges)))
    pool.close()

def col(edges):
    pool = multiprocessing.Pool(1)
    pool.map(color, range(0, len(edges)))
    pool.close()
    subprocess.check_output('samtools merge %s/bam/coloredBAM.bam -f `find %s/bam -name "*unitig*.bam"`' % (output,output), shell=True, capture_output=False)
    subprocess.check_output('rm `find %s/bam -name "*unitig*.bam"`' % output, shell=True,capture_output=False)
    pysam.samtools.index("%s/bam/coloredBAM.bam" % output, "%s/bam/coloredBAM.bai" % output)

import numpy as np
#all_data={}
#np.save("output/all_data.npy", all_data)


if __name__ == "__main__":
    phase(edges)
    col(edges)

#np.save("output/all_data.npy", all_data)




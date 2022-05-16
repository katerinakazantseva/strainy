import multiprocessing
from params import *
from cluster import cluster
from color_bam import color
import pysam
import os
import subprocess


dirs = ("output/vcf/","output/adj_M/","output/clusters/","output/graphs/","output/bam/")
for dir in dirs:
    try:
        os.stat(dir)
    except:
         os.mkdir(dir)

def phase(edges):
    pool = multiprocessing.Pool(3)
    pool.map(cluster, range(0, len(edges)))
    pool.close()

def col(edges):
    pool = multiprocessing.Pool(3)
    pool.map(color, range(0, len(edges)))
    pool.close()
    subprocess.check_output('samtools merge output/bam/coloredBAM.bam -f `find output/bam -name "*edge*.bam"`', shell=True, capture_output=False)
    subprocess.check_output('rm `find output/bam -name "*edge*.bam"`', shell=True,capture_output=False)
    pysam.samtools.index("output/bam/coloredBAM.bam", "output/bam/coloredBAM.bai")




if __name__ == "__main__":
    phase(edges)
    #col(edges)
    #transform(edges)



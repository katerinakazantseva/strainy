import multiprocessing
import pysam
import pickle
import os
import subprocess
from multiprocessing.managers import BaseManager

from cluster import cluster
from color_bam import color
from flye_consensus import FlyeConsensus
from params import *





dirs = ( "%s/vcf/" % output ,"%s/adj_M/" % output,"%s/clusters/" % output,"%s/graphs/" % output,"%s/bam/" % output,"%s/bam/clusters" % output,   "%s/flye_inputs" % output,"%s/flye_outputs" % output)

for dir in dirs:
    try:
        os.stat(dir)
    except:

        os.makedirs(dir)



class CustomManager(BaseManager):
    pass


def phase(edges):
    CustomManager.register('FlyeConsensus', FlyeConsensus)

    with CustomManager() as manager:
        default_manager = multiprocessing.Manager()
        lock = default_manager.Lock()
        empty_consensus_dict = {}
        num_processes = multiprocessing.cpu_count() if processes == -1 else processes
        shared_flye_consensus = manager.FlyeConsensus(bam, gfa, num_processes, empty_consensus_dict, lock)
        pool = multiprocessing.Pool(num_processes)
        init_args = [(i, shared_flye_consensus) for i in range(len(edges))]
        pool.map(cluster, init_args)
        pool.close()
        shared_flye_consensus.print_cache_statistics()
        return shared_flye_consensus.get_consensus_dict()




def col(edges):
    pool = multiprocessing.Pool(1)
    pool.map(color, range(0, len(edges)))
    pool.close()
    subprocess.check_output('samtools merge %s/bam/coloredBAM.bam -f `find %s/bam -name "*unitig*.bam"`' % (output,output), shell=True, capture_output=False)
    subprocess.check_output('rm `find %s/bam -name "*unitig*.bam"`' % output, shell=True,capture_output=False)
    pysam.samtools.index("%s/bam/coloredBAM.bam" % output, "%s/bam/coloredBAM.bai" % output)

import numpy as np


if __name__ == "__main__":
    consensus_dict = phase(edges)
    if write_consensus_cache:
        with open(consensus_cache_path, 'wb') as f:
            pickle.dump(consensus_dict, f)
    #col(edges)

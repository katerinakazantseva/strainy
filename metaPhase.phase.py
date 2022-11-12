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



dirs = ("output/",
        "output/vcf/",
        "output/adj_M/",
        "output/clusters/",
        "output/graphs/",
        "output/bam/",
        "output/bam/clusters/",
        "output/flye_inputs",
        "output/flye_outputs")


for dir in dirs:
    try:
        os.stat(dir)
    except:
        os.mkdir(dir)



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
    #subprocess.check_output('samtools merge output/bam/coloredBAM.bam -f `find output/bam -name "*edge*.bam"`',
                            shell=True, capture_output=False)
    #subprocess.check_output('rm `find output/bam -name "*edge*.bam"`', shell=True, capture_output=False)
    #pysam.index("output/bam/coloredBAM.bam")


if __name__ == "__main__":
    consensus_dict = phase(edges)
    if write_consensus_cache:
        with open(consensus_cache_path, 'wb') as f:
            pickle.dump(consensus_dict, f)
    col(edges)

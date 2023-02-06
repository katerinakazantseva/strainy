import multiprocessing
import pysam
import pickle
import os
import sys
import subprocess
from multiprocessing.managers import BaseManager
import logging
import shutil

from strainy.clustering.cluster import cluster
from strainy.color_bam import color
from strainy.flye_consensus import FlyeConsensus
from strainy.params import *
from strainy.logging import set_thread_logging


logger = logging.getLogger()



def _thread_fun(args):
    set_thread_logging(StRainyArgs.log_phase, "phase", multiprocessing.current_process().pid)
    logger.info("\n\n\t==== Processing uniting " + str(StRainyArgs.edges[args[0]]) + " ====")
    cluster(*args)
    logger.info("Thread finished!")


def _error_callback(pool, e):
    logger.error("Worker thread exception! " + str(e))
    pool.terminate()
    raise e


def phase(edges):
    logger.info("CMD: " + " ".join(sys.argv[1:]))
    if os.path.isdir(StRainyArgs.log_phase):
        shutil.rmtree(StRainyArgs.log_phase)
    os.mkdir(StRainyArgs.log_phase)
    #CustomManager.register('FlyeConsensus', FlyeConsensus)
    #with CustomManager() as manager:
    default_manager = multiprocessing.Manager()
    #lock = default_manager.Lock()
    empty_consensus_dict = {}
    num_processes = multiprocessing.cpu_count() if StRainyArgs.threads == -1 else StRainyArgs.threads
    #shared_flye_consensus = manager.FlyeConsensus(StRainyArgs.bam, StRainyArgs.gfa, num_processes, empty_consensus_dict, lock)
    shared_flye_consensus = FlyeConsensus(StRainyArgs.bam, StRainyArgs.fa, num_processes, empty_consensus_dict, default_manager)
    pool = multiprocessing.Pool(num_processes)
    init_args = [(i, shared_flye_consensus) for i in range(len(edges))]
    pool.map_async(_thread_fun, init_args, error_callback=lambda e: _error_callback(pool, e))
    pool.close()
    pool.join()
    shared_flye_consensus.print_cache_statistics()
    return shared_flye_consensus.get_consensus_dict()


def color_bam(edges):
    pool = multiprocessing.Pool(3)
    pool.map(color, range(0, len(edges)))
    pool.close()

    out_bam_dir = os.path.join(StRainyArgs.output, 'bam')
    to_merge_file = os.path.join(out_bam_dir, 'to_merge.txt')
    to_delete = []
    with open(to_merge_file, "wb") as f:
        for fname in subprocess.check_output('find %s -name "*unitig*.bam"' % out_bam_dir, shell=True).split(b'\n'):
            if len(fname):
                f.write(fname + b'\n')
                to_delete.append(fname)

    subprocess.check_output('samtools merge %s/bam/coloredBAM.bam -f -b %s' % (StRainyArgs.output, to_merge_file),
                            shell=True, capture_output=False)
    for f in to_delete:
        os.remove(f)
    #subprocess.check_output('rm `find %s/bam -name "*unitig*.bam"`' % StRainyArgs.output, shell=True, capture_output=False)
    pysam.samtools.index("%s/bam/coloredBAM.bam" % StRainyArgs.output, "%s/bam/coloredBAM.bai" % StRainyArgs.output)


def phase_main():
    #logging.info(StRainyArgs)
    logger.info("Starting phasing")
    dirs = ("%s/vcf/" % StRainyArgs.output,
            "%s/adj_M/" % StRainyArgs.output,
            "%s/clusters/" % StRainyArgs.output,
            "%s/graphs/" % StRainyArgs.output,
            "%s/bam/" % StRainyArgs.output,
            "%s/bam/clusters" % StRainyArgs.output,
            "%s/flye_inputs" % StRainyArgs.output,
            "%s/flye_outputs" % StRainyArgs.output)

    for dir in dirs:
        try:
            os.stat(dir)
        except:
            os.makedirs(dir)

    consensus_dict = phase(StRainyArgs.edges)
    if write_consensus_cache:
        with open(os.path.join(StRainyArgs.output, consensus_cache_path), 'wb') as f:
            pickle.dump(consensus_dict, f)
    color_bam(StRainyArgs.edges)


if __name__ == "__main__":
    phase_main()

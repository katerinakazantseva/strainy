import multiprocessing
import pysam
import pickle
import os
import sys
import subprocess
import multiprocessing
import logging
import shutil
import traceback
import time

from strainy.clustering.cluster import cluster
from strainy.color_bam import color
from strainy.flye_consensus import FlyeConsensus
from strainy.params import *
from strainy.logging import set_thread_logging


logger = logging.getLogger()


def _thread_fun(i, shared_flye_consensus, args):
    init_global_args_storage(args)

    set_thread_logging(StRainyArgs().log_phase, "phase", multiprocessing.current_process().pid)
    logger.info("\n\n\t == == Processing unitig " + str(StRainyArgs().edges_to_phase[i]) + " == == ")

    try:
        cluster(i, shared_flye_consensus)
    except Exception as e:
        logger.error("Worker thread exception! " + str(e) + "\n" + traceback.format_exc())
        raise e

    logger.debug("Thread worker function finished!")


def phase(edges, args):
    logger.info("CMD: " + " ".join(sys.argv[1:]))

    empty_consensus_dict = {}
    default_manager = multiprocessing.Manager()
    shared_flye_consensus = FlyeConsensus(StRainyArgs().bam, StRainyArgs().fa, 1, empty_consensus_dict, default_manager)
    if StRainyArgs().threads == 1:
        for i in range(len(edges)):
            cluster(i, shared_flye_consensus)
    else:
        pool = multiprocessing.Pool(StRainyArgs().threads)
        init_args = [(i, shared_flye_consensus, args) for i in range(len(edges))]

        results = pool.starmap_async(_thread_fun, init_args, chunksize=1)
        while not results.ready():
            time.sleep(0.01)
            if not results._success:
                pool.terminate()
                raise Exception("Error in worker thread, exiting")

        pool.close()
        pool.join()

    shared_flye_consensus.print_cache_statistics()
    return shared_flye_consensus.get_consensus_dict()


def color_bam(edges):
    logger.info("Creating phased bam")
    for e in edges:
        color(e)

    out_bam_dir = os.path.join(StRainyArgs().output_intermediate, "bam")
    final_aln = os.path.join(StRainyArgs().output, "alignment_phased.bam")

    files_to_be_merged = []
    for fname in subprocess.check_output(f'find {out_bam_dir} -name "*unitig*.bam"', shell = True, universal_newlines = True).split("\n"):
        if len(fname):
            files_to_be_merged.append(fname)

    # Number of file to be merged could be > 4092, in which case samtools merge throws too many open files error
    for i, bam_file in enumerate(files_to_be_merged):
        # fetch the header and put it at the top of the file, for the first bam_file only
        if i == 0:
            subprocess.check_output(f'samtools view -H {bam_file} > {out_bam_dir}/coloredSAM.sam', shell = True)

        # convert bam to sam, append to the file
        subprocess.check_output(f'samtools view {bam_file} >> {out_bam_dir}/coloredSAM.sam', shell = True)

    # convert the file to bam and sort
    subprocess.check_output(f'samtools view -b {out_bam_dir}/coloredSAM.sam >> {out_bam_dir}/unsortedBAM.bam', shell = True)
    pysam.samtools.sort(f'{out_bam_dir}/unsortedBAM.bam', "-o", final_aln)
    pysam.samtools.index(final_aln)

    # remove unnecessary files
    os.remove(f'{out_bam_dir}/unsortedBAM.bam')
    os.remove(f'{out_bam_dir}/coloredSAM.sam')
    for f in files_to_be_merged:
        os.remove(f)


def phase_main(args):
    #logging.info(StRainyArgs)
    logger.info("Starting phasing")
    dirs = ("%s/vcf/" % StRainyArgs().output_intermediate,
            "%s/clusters/" % StRainyArgs().output_intermediate,
            "%s/bam/" % StRainyArgs().output_intermediate,
            "%s/bam/clusters" % StRainyArgs().output_intermediate,
            "%s/flye_inputs" % StRainyArgs().output_intermediate,
            "%s/flye_outputs" % StRainyArgs().output_intermediate
    )
    
    debug_dirs = ("%s/graphs/" % StRainyArgs().output_intermediate,
                  "%s/adj_M/" % StRainyArgs().output_intermediate
    )

    for dir in dirs:
        os.makedirs(dir, exist_ok=True)

    if StRainyArgs().debug:
        os.makedirs(dir, exist_ok=True)

    consensus_dict = phase(StRainyArgs().edges_to_phase, args)
    if write_consensus_cache:
        with open(os.path.join(StRainyArgs().output_intermediate, consensus_cache_path), "wb") as f:
            pickle.dump(consensus_dict, f)
    color_bam(StRainyArgs().edges)
    logger.info("Done")


if __name__ == "__main__":
    phase_main()

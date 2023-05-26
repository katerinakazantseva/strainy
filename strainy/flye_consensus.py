import gzip
import multiprocessing
import subprocess
import os
import shutil
import random
import logging
import time

import edlib
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import Namespace

from strainy.params import *
#from flye.main import _run_polisher_only

logger = logging.getLogger()
logging.basicConfig(level=logging.DEBUG)

def calculate_coverage(position, bed_file_content):
    """
    Calculates and returns the coverage for a given position that is relative to the reference seq, not the aligment
    string
    """
    for row in bed_file_content:
        # row = [interval_start, interval_end, coverage]
        if row[0] <= position <= row[1]:
            return row[2]
    logger.warning("Coordinate not found in .bed file, assuming coverage is 0")
    return 0


class FlyeConsensus:
    def __init__(self, bam_file_name, graph_fasta_name, num_processes, consensus_dict, multiproc_manager,
                indel_block_length_leniency=5):

        self._consensus_dict = multiproc_manager.dict(consensus_dict)
        self._lock = multiproc_manager.Lock()

        self._bam_path = bam_file_name
        self._read_index = None
        self._unitig_seqs = {}
        for seq in SeqIO.parse(graph_fasta_name, "fasta"):
            self._unitig_seqs[str(seq.id)] = str(seq.seq)

        self._num_processes = num_processes
        self._indel_block_length_leniency = indel_block_length_leniency
        self._coverage_limit = min_consensus_cov[StRainyArgs().mode]
        if StRainyArgs().mode == "hifi":
            self._platform = "pacbio"
            self._read_type = "hifi"
            self._mode = "--pacbio-hifi"
        elif StRainyArgs().mode == "nano":
            self._platform = "nano"
            self._read_type = "raw"
            self._mode = "--nano-raw"

        self._key_hit = multiproc_manager.Value("i", 0)
        self._key_miss = multiproc_manager.Value("i", 0)

        self._debug_count = multiproc_manager.Value("i", 0)
        self._call_count = multiproc_manager.Value("i", 0)

    def get_consensus_dict(self):
        return self._consensus_dict.copy()

    def print_cache_statistics(self):
        logger.info(f"Total number of key hits and misses for consensus computation:")
        logger.info(f" H:{self._key_hit.value}, M:{self._key_miss.value}")


    def extract_reads(self, read_names, output_file, edge=""):
        """
        based on the code by Tim Stuart https://timoast.github.io/blog/2015-10-12-extractreads/
        Extract the reads given query names to a new bam file
        """
        cluster_start = -1
        cluster_end = -1
        read_limits = []

        read_list = []  # stores the reads to be written after the cluster start/end is calculated

        #ts = time.time()
        if self._read_index is None:
            self._read_index = pysam.IndexedReads(pysam.AlignmentFile(self._bam_path, "rb"))
            self._read_index.build()
        #te = time.time()
        #logger.debug("Index building time %f", te - ts)

        for name in read_names:
            iterator = self._read_index.find(name)
            for x in iterator:
                if x.reference_name == edge:
                    if x.reference_start < cluster_start or cluster_start == -1:
                        cluster_start = x.reference_start
                    if x.reference_end > cluster_end or cluster_end == -1:
                        cluster_end = x.reference_end
                    read_list.append(x)
                    read_limits.append((x.reference_start, x.reference_end))

        #out = pysam.Samfile(output_file, "wb", header=self._bam_header)
        out = pysam.Samfile(output_file, "wb", template=pysam.AlignmentFile(self._bam_path, "rb"))
        for x in read_list:
            temp_dict = x.to_dict()
            temp_dict["ref_pos"] = str(int(temp_dict["ref_pos"]) - cluster_start)
            y = x.from_dict(temp_dict, x.header)  # create a new read from the modified dictionary
            out.write(y)

        out.close()

        return cluster_start, cluster_end, read_limits

    def flye_consensus(self, cluster, edge, cl, debug=False):
        """
        Computes the Flye based consensus of a cluster of reads for a specific edge.
        cluster: id (int)
        cl: dataframe with columns read_name and cluster(id)
        edge: edge name (str)
        """

        # check if the output for this cluster-edge pair exists in the cache
        consensus_dict_key = f"{cluster}-{edge}"
        with self._lock:
            if consensus_dict_key in self._consensus_dict:
                self._key_hit.value += 1
                return self._consensus_dict[consensus_dict_key]
            self._key_miss.value += 1

        # fetch the read names in this cluster and extract those reads to a new bam file to be used by the
        # Flye polisher
        reads_from_curr_cluster = cl.loc[cl["Cluster"] == cluster]["ReadName"].to_numpy()  # store read names
        salt = random.randint(1000, 10000)
        fprefix = "%s/flye_inputs/" % StRainyArgs().output
        cluster_start, cluster_end, read_limits = self.extract_reads(reads_from_curr_cluster,
                                                        f"{fprefix}cluster_{cluster}_reads_{salt}.bam", edge)

        logger.debug((f"CLUSTER:{cluster}, CLUSTER_START:{cluster_start}, CLUSTER_END:{cluster_end}, EDGE:{edge},"
               f"# OF READS:{len(reads_from_curr_cluster)}"))

        # access the edge in the graph and cut its sequence according to the cluster start and end positions
        # this sequence is written to a fasta file to be used by the Flye polisher
        edge_seq_cut = self._unitig_seqs[edge][cluster_start:cluster_end]
        fname = f"{fprefix}{edge}-cluster{cluster}-{salt}"
        record = SeqRecord(
            Seq(edge_seq_cut),
            id=f"{edge}",
            name=f"{edge} sequence cut for cluster {cluster}",
            description=""
        )
        SeqIO.write([record], f"{fname}.fa", "fasta")

        # sort the bam file
        pysam.sort("-o", f"{fprefix}cluster_{cluster}_reads_sorted_{salt}.bam",
                   f"{fprefix}cluster_{cluster}_reads_{salt}.bam")
        # index the bam file
        pysam.index(f"{fprefix}cluster_{cluster}_reads_sorted_{salt}.bam")

        #  Polisher arguments for to call _run_polisher_only(polish_args)
        # polish_args = Namespace(polish_target=f"{fname}.fa",
        #                         reads=[f"{fprefix}cluster_{cluster}_reads_sorted_{salt}.bam"],
        #                         out_dir=f"{MetaPhaseArgs.output}/flye_outputs/flye_consensus_{edge}_{cluster}_{salt}",
        #                         num_iters=1,
        #                         threads=1,
        #                         platform=self._platform,
        #                         read_type=self._read_type,
        #                         output_progress=True)

        polish_cmd = f"{StRainyArgs().flye} --polish-target {fname}.fa --threads {self._num_processes}" \
                     f" {self._mode} {fprefix}cluster_{cluster}_reads_sorted_{salt}.bam " \
                     f"-o {StRainyArgs().output}/flye_outputs/flye_consensus_{edge}_{cluster}_{salt}"
        try:
            logger.debug("Running Flye polisher")
            subprocess.check_output(polish_cmd, shell=True, capture_output=False, stderr=open(os.devnull, "w"))
            logger.debug("Running Flye polisher - finished!")
        except subprocess.CalledProcessError as e:
            logger.error("Error running the Flye polisher. Make sure the fasta file contains only the primary alignments")
            logger.error(e)
            with self._lock:
                self._consensus_dict[consensus_dict_key] = {
                    'consensus': Seq(''),
                    'start': cluster_start,
                    'end': cluster_end
                }
                return self._consensus_dict[consensus_dict_key]

        try:
            # read back the output of the Flye polisher
            consensus = SeqIO.read(f"{StRainyArgs().output}/flye_outputs/flye_consensus_{edge}_{cluster}_{salt}/polished_1.fasta",
                                   "fasta")
        except (ImportError, ValueError) as e:
            # If there is an error, the sequence string is set to empty by default
            logger.warning("WARNING: error reading back the flye output, defaulting to empty sequence for consensus")
            if type(e).__name__ == 'ImportError':
                logger.warning('found ImportError')
            consensus = SeqRecord(
                seq=''
            )
        # delete the created input files to Flye
        if delete_files:
            os.remove(f"{fname}.fa")
            os.remove(f"{fprefix}cluster_{cluster}_reads_{salt}.bam")
            os.remove(f"{fprefix}cluster_{cluster}_reads_sorted_{salt}.bam")
            os.remove(f"{fprefix}cluster_{cluster}_reads_sorted_{salt}.bam.bai")
            shutil.rmtree(f"{StRainyArgs().output}/flye_outputs/flye_consensus_{edge}_{cluster}_{salt}")

        with self._lock:
            self._consensus_dict[consensus_dict_key] = {
                'consensus': consensus.seq,
                'start': cluster_start,
                'end': cluster_end,
                'read_limits': read_limits,
                'bam_path': f"{fprefix}cluster_{cluster}_reads_{salt}.bam",
                'reference_path': f"{fname}.fa",
                'bed_path': f"{StRainyArgs().output}/flye_outputs/flye_consensus_{edge}_{cluster}_{salt}/"
                            f"base_coverage.bed.gz"
            }
            return self._consensus_dict[consensus_dict_key]

    def _edlib_align(self, seq_a, seq_b):
        band_size = 32
        aln = None
        while True:
            aln = edlib.align(seq_a, seq_b, "NW", "path", band_size)
            if aln["editDistance"] == -1:   #need to increase alignment band
                if band_size > max(len(seq_a), len(seq_b)):
                    raise Exception("Something's wrong with edlib")
                band_size *= 2
            else:
                break
        #print("Band size", band_size)
        nice = edlib.getNiceAlignment(aln, seq_a, seq_b)
        #note that target and query are swapped because the definition is edlib.align(query, target)
        return nice["query_aligned"], nice["target_aligned"], nice["matched_aligned"]

    def _parse_bed_coverage(self, filename):
        contents = []
        with gzip.open(filename, 'r') as f:
            for line in f:
                try:
                    contents.append(list(map(int, line.strip().split()[1:])))
                except ValueError:
                    pass
        return contents

    def _custom_scoring_function(self, aligned_first, alignment_string, aligned_second, intersection_start, first_cl_dict, second_cl_dict):
        """
        A custom distance scoring function for two sequences taking into account the artifacts of Flye consensus.
        alignment_string: a string consisting of '-', '.', '|' characters which correspond to indel, mismatch, match,
        respectively.
        Mismatches are worth 1 point, indels are 1 point each if there are more than 5 of them in a contiguous block.
        Except for the indel blocks start at the beginning or finish at the end, those indels are ignored.
        Moreover, variants that are covered by less than self._coverage_limits are ignored (assumed match)
        """
        score = 0
        indel_length = 0
        indel_block_start = -1
        alignment_list = list(alignment_string)

        first_shift = intersection_start - first_cl_dict['start']
        second_shift = intersection_start - second_cl_dict['start']

        # read the contents of the bed.gz files that contain coverage information for each coordinate
        # list containing lists with (start, end, coverage) for each coordinate interval
        cl1_bed_contents = self._parse_bed_coverage(first_cl_dict['bed_path'])
        cl2_bed_contents = self._parse_bed_coverage(second_cl_dict['bed_path'])

        for i in range(len(alignment_list)):
            if alignment_list[i] not in "-.|":
                raise Exception("Unknown alignment sybmol!")

            # true coordinate = current coordinate on the alignment_string
            # + start of the intersection
            # - gaps in the target (or query) sequence thus far
            cl1_true_coor = first_shift + i - aligned_first[:i].count('-')
            cl2_true_coor = second_shift + i - aligned_second[:i].count('-')

            # ignore variants with low coverage
            if ((alignment_list[i] == '-' or alignment_list[i] == '.')
                    and (calculate_coverage(cl1_true_coor, cl1_bed_contents) < self._coverage_limit
                    or calculate_coverage(cl2_true_coor, cl2_bed_contents) < self._coverage_limit)):
                alignment_list[i] = '|'

            if alignment_list[i] == '-':
                # igonre the indels at the first or last position
                if i == 0:
                    indel_length = 1
                    indel_block_start = i
                    continue
                elif i == (len(alignment_list) - 1):
                    continue
                else:
                    # start of an indel block
                    if alignment_list[i - 1] != '-':
                        # an indel blocks starting at position 0 will be ignored
                        indel_block_start = i
                        indel_length = 1
                    # continuing indel block
                    else:
                        indel_length += 1

            # a contiguous gap ends
            elif alignment_list[i - 1] == '-':
                if indel_length >= self._indel_block_length_leniency and indel_block_start != 0:
                    score += indel_length
                indel_length = 0

            # mismatch
            if alignment_list[i] == '.':
                score += 1
        return score

    def _log_alignment_info(self, aligned_first, alignment_string, aligned_second, first_cl_dict, second_cl_dict,
                            score, intersection_start, intersection_end):
        #if self._num_processes == 1:
        #    pid = 0
        #else:
        #    pid = int(re.split('\'|-', str(multiprocessing.current_process()))[2]) - 1
        pid = int(multiprocessing.current_process().pid)
        fname = f"{StRainyArgs().output}/distance_inconsistency-{pid}.log"

        first_shift = intersection_start - first_cl_dict['start']
        second_shift = intersection_start - second_cl_dict['start']
        cl1_bed_contents = self._parse_bed_coverage(first_cl_dict['bed_path'])
        cl2_bed_contents = self._parse_bed_coverage(second_cl_dict['bed_path'])

        with open(fname, 'a+') as f:
            """
            Things to write to the file:
            A different file for each process (filename{pid}.log)
            Alignment string,
            calculated score
            for each cluster:
                consensuses
                read start and end positions 
            """
            f.write("ALIGNMENT:\n")
            f.write(alignment_string + '\n')
            mismatch_positions = [i for i in range(len(alignment_string)) if alignment_string.startswith('.', i)]
            ref_mismatch = [i + intersection_start for i in mismatch_positions]
            f.write(f"# OF MISMATCHES: {len(mismatch_positions)}\n")
            f.write(f"MISMATCH POSITIONS: {mismatch_positions}\n")
            f.write(f"MISMATCH POSITIONS IN REF (APPROXIMATE): {ref_mismatch}\n")
            f.write(f"SCORE:{score}\n")
            length = intersection_end - intersection_start
            f.write(f"INTERSECTION AREA:({intersection_start}, {intersection_end}),"
                    f"LENGTH:{length}\n")

            # Calculate coverages of mismatches for both clusters
            first_cl_mismatch_coverages = []
            second_cl_mismatch_coverages = []
            for i in mismatch_positions:
                # i is relative to the alignment string
                cl1_true_coor = first_shift + i - aligned_first[:i].count('-')
                cl2_true_coor = second_shift + i - aligned_second[:i].count('-')

                first_cl_mismatch_coverages.append(
                    calculate_coverage(cl1_true_coor, cl1_bed_contents)
                )
                second_cl_mismatch_coverages.append(
                    calculate_coverage(cl2_true_coor, cl2_bed_contents)
                )

            f.write("FIRST CLUSTER:\n")
            f.write("\tREADS:\n")
            f.write(f"\t{first_cl_dict['read_limits']}\n")
            f.write(f"\tMISMATCH COVERAGES: {first_cl_mismatch_coverages}\n")
            avg_coverage = 0
            for read in first_cl_dict['read_limits']:
                avg_coverage += read[1] - read[0]
            avg_coverage = "{:.2f}".format(avg_coverage / len(first_cl_dict['consensus']))
            f.write(f"\tAVERAGE COVERAGE:{avg_coverage}\n")
            f.write(f"\tBAM FILE PATH:{first_cl_dict['bam_path']}\n")
            f.write(f"\tREFERENCE FILE PATH:{first_cl_dict['reference_path']}\n")

            f.write("SECOND CLUSTER:\n")
            f.write("\tREADS:\n")
            f.write(f"\t{second_cl_dict['read_limits']}\n")
            f.write(f"\tMISMATCH COVERAGES: {second_cl_mismatch_coverages}\n")
            avg_coverage = 0
            for read in second_cl_dict['read_limits']:
                avg_coverage += read[1] - read[0]
            avg_coverage = "{:.2f}".format(avg_coverage / len(second_cl_dict['consensus']))
            f.write(f"\tAVERAGE COVERAGE:{avg_coverage}\n")
            f.write(f"\tBAM FILE PATH:{second_cl_dict['bam_path']}\n")
            f.write(f"\tREFERENCE FILE PATH:{second_cl_dict['reference_path']}\n")

            f.write("**********-------************\n\n")


    def cluster_distance_via_alignment(self, first_cl, second_cl, cl, edge, debug=False):
        """
        Computes the distance between two clusters consensus'. The distance is based on the global alignment between the
        intersecting parts of the consensus'.
        first_cl: id (int)
        second_cl: id (int)
        cl: dataframe with columns read_name and cluster(id)
        edge: edge name (str)
        """
        self._call_count.value += 1
        if debug:
            self._debug_count.value += 1
        if self._debug_count.value > 0:
            logger.debug(f"{self._debug_count.value}/{self._call_count.value} disagreements")
        first_cl_dict = self.flye_consensus(first_cl, edge, cl, debug)
        second_cl_dict = self.flye_consensus(second_cl, edge, cl, debug)

        intersection_start = max(first_cl_dict['start'], second_cl_dict['start'])
        intersection_end = min(first_cl_dict['end'], second_cl_dict['end'])

        # clip the intersecting parts of both consensus'
        first_consensus_clipped = first_cl_dict['consensus'][
                                  intersection_start - first_cl_dict['start']:intersection_end - first_cl_dict['start']]
        second_consensus_clipped = second_cl_dict['consensus'][
                                   intersection_start - second_cl_dict['start']:intersection_end - second_cl_dict['start']]

        if (intersection_end - intersection_start < 1
                or len(first_consensus_clipped) == 0
                or len(second_consensus_clipped) == 0):
            logger.debug(f'Intersection length for clusters is less than 1 for clusters {first_cl}, {second_cl} in {edge}')
            return 1


        aligned_first, aligned_second, edlib_aln = self._edlib_align(first_consensus_clipped, second_consensus_clipped)
        edlib_score = self._custom_scoring_function(aligned_first, edlib_aln, aligned_second, intersection_start,
                                                    first_cl_dict, second_cl_dict)

        if debug:
            self._log_alignment_info(aligned_first, edlib_aln, aligned_second,
                                    first_cl_dict, second_cl_dict,
                                    edlib_score,
                                    intersection_start, intersection_end)

        # score is not normalized!
        return edlib_score

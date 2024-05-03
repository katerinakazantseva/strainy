import gzip
import traceback
import subprocess
import os
import shutil
import random
import logging
import sys

import edlib
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import Namespace

from strainy.params import *

logger = logging.getLogger()
logging.basicConfig(level=logging.DEBUG)

from flye.main import _run_polisher_only

def calculate_coverage(position, bed_file_content):
    """
    Calculates and returns the coverage for a given position that is relative to the reference seq, not the aligment
    string
    """
    for row in bed_file_content:
        # row = [interval_start, interval_end, coverage]
        if row[0] <= position <= row[1]:
            return row[2]
    logger.debug("Coordinate not found in .bed file, assuming coverage is 0")
    return 0


class FlyeConsensus:
    def __init__(self, bam_file_name, graph_fasta_name, num_processes, consensus_dict, multiproc_manager,
                indel_block_length_leniency=5):

        self._lock = multiproc_manager.Lock()

        self._consensus_dict = multiproc_manager.dict(consensus_dict)
        self._alignment_cache= multiproc_manager.dict()

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


        self._position_hit = multiproc_manager.Value("i", 0)
        self._position_miss = multiproc_manager.Value("i", 0)

        self._alignment_cache_hit = multiproc_manager.Value("i", 0)
        self._alignment_cache_miss = multiproc_manager.Value("i", 0)


    def get_consensus_dict(self):
        return self._consensus_dict.copy()


    def print_cache_statistics(self):
        logger.info(f"Total number of key hits and misses for consensus computation:")
        logger.info(f" H:{self._key_hit.value}, M:{self._key_miss.value}")
        logger.info(f"Position hit/miss")
        logger.info(f" H:{self._position_hit.value}, M:{self._position_miss.value}")
        logger.info(f"Alignment cache hit/miss")
        logger.info(f" H:{self._alignment_cache_hit.value}, M:{self._alignment_cache_miss.value}")


    def _extract_reads(self, read_names, start_pos, output_file, edge=""):
        """
        based on the code by Tim Stuart https://timoast.github.io/blog/2015-10-12-extractreads/
        Extract the reads given query names to a new bam file
        """
        cluster_start = -1
        cluster_end = -1
        read_limits = []

        read_list = []  # stores the reads to be written after the cluster start/end is calculated

        if self._read_index is None:
            self._read_index = pysam.IndexedReads(pysam.AlignmentFile(self._bam_path, "rb"))
            self._read_index.build()

        for i, name in enumerate(read_names):
            iterator = self._read_index.find(name)
            for x in iterator:
                if x.reference_name == edge and x.reference_start == start_pos[i]:
                    if x.reference_start < cluster_start or cluster_start == -1:
                        cluster_start = x.reference_start
                    if x.reference_end > cluster_end or cluster_end == -1:
                        cluster_end = x.reference_end
                    read_list.append(x)
                    read_limits.append((x.reference_start, x.reference_end))

        out = pysam.Samfile(output_file, "wb", template=pysam.AlignmentFile(self._bam_path, "rb"))
        for x in read_list:
            temp_dict = x.to_dict()
            temp_dict["ref_pos"] = str(int(temp_dict["ref_pos"]) - cluster_start)
            y = x.from_dict(temp_dict, x.header)  # create a new read from the modified dictionary
            out.write(y)

        out.close()

        return cluster_start, cluster_end, read_limits
    

    def _clip_consensus_seq(self, sequence, read_limits, bed_contents, curr_start, coverage_limit):
        start_pos = sorted([start for (start, _) in read_limits])
        end_pos = sorted([end for (_, end) in read_limits], reverse=True)

        try:
            new_start_pos = start_pos[coverage_limit - 1]
            new_end_pos = end_pos[coverage_limit - 1]
        except IndexError:
            new_start_pos = start_pos[0]
            new_end_pos = end_pos[0]

        return new_start_pos, new_end_pos, sequence[new_start_pos - curr_start:new_end_pos - curr_start]


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
        start_pos_of_reads = cl.loc[cl["Cluster"] == cluster]["Start"].to_numpy()
        salt = random.randint(1000, 10000)
        fprefix = "%s/flye_inputs/" % StRainyArgs().output_intermediate
        bam_subset = f"{fprefix}{edge}_cluster_{cluster}_reads_{salt}.bam"
        bam_subset_sorted = f"{fprefix}{edge}_cluster_{cluster}_reads_{salt}_sorted.bam"
        cluster_start, cluster_end, read_limits = self._extract_reads(reads_from_curr_cluster, start_pos_of_reads,
                                                                      bam_subset, edge)

        logger.debug((f"CLUSTER:{cluster}, CLUSTER_START:{cluster_start}, CLUSTER_END:{cluster_end}, EDGE:{edge},"
               f"# OF READS:{len(reads_from_curr_cluster)}"))

        # access the edge in the graph and cut its sequence according to the cluster start and end positions
        # this sequence is written to a fasta file to be used by the Flye polisher
        ref_seq_cut = self._unitig_seqs[edge][cluster_start:cluster_end]
        fname = f"{fprefix}{edge}-cluster{cluster}-{salt}"
        record = SeqRecord(
            Seq(ref_seq_cut),
            id=f"{edge}",
            name=f"{edge} sequence cut for cluster {cluster}",
            description=""
        )
        SeqIO.write([record], f"{fname}.fa", "fasta")

        try:
            # sort the bam file
            pysam.sort("-o", bam_subset_sorted, bam_subset)
            # index the bam file
            pysam.index(bam_subset_sorted)
        except pysam.utils.SamtoolsError  as e:
            logger.error(f'Error while sorting {bam_subset_sorted}')
            logger.error(traceback.format_exc())

        #  Polisher arguments for to call _run_polisher_only(polish_args)
        flye_out_dir = f"{StRainyArgs().output_intermediate}/flye_outputs/flye_consensus_{edge}_{cluster}_{salt}"
        polish_args = Namespace(polish_target=f"{fname}.fa",
                                reads=[bam_subset_sorted],
                                out_dir=flye_out_dir,
                                num_iters=1,
                                threads=1,
                                platform=self._platform,
                                read_type=self._read_type)

        # polish_cmd = f"{StRainyArgs().flye} --polish-target {fname}.fa --threads {self._num_processes}" \
        #              f" {self._mode} {fprefix}cluster_{cluster}_reads_sorted_{salt}.bam " \
        #              f"-o {StRainyArgs().output_intermediate}/flye_outputs/flye_consensus_{edge}_{cluster}_{salt}"
        try:
            logger.debug("Running Flye polisher")
            # subprocess.check_output(polish_cmd, shell=True, capture_output=False, stderr=open(os.devnull, "w"))
            # TODO: this should move to the top when flye pull request is merged
            if not os.path.isdir(polish_args.out_dir):
                os.mkdir(polish_args.out_dir)
            _run_polisher_only(polish_args, output_progress=False)
            logger.debug("Running Flye polisher - finished!")
        except Exception as e:
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
            consensus = SeqIO.read(os.path.join(flye_out_dir, "polished_1.fasta"), "fasta")
        except (ImportError, ValueError) as e:
            # If there is an error, the sequence string is set to empty by default
            logger.warning("WARNING: error reading back the flye output, defaulting to empty sequence for consensus")
            if type(e).__name__ == 'ImportError':
                logger.warning('found ImportError')
            consensus = SeqRecord(
                seq=''
            )

        bed_content = self._parse_bed_coverage(os.path.join(flye_out_dir, "base_coverage.bed.gz"))

        # delete the created input files to Flye
        if delete_flye_files:
            try:
                os.remove(f"{fname}.fa")
                os.remove(bam_subset)
                os.remove(bam_subset_sorted)
                os.remove(bam_subset_sorted + ".bai")
                shutil.rmtree(flye_out_dir)
            except (OSError, FileNotFoundError):
                pass

        start, end, consensus_clipped = self._clip_consensus_seq(consensus.seq, 
                                                                 read_limits, 
                                                                 bed_content,
                                                                 cluster_start,
                                                                 2)
        with self._lock:
            self._consensus_dict[consensus_dict_key] = {
                'consensus': consensus_clipped,
                'start': start,
                'end': end,
                'read_limits': read_limits,
                'bam_path': bam_subset,
                'reference_path': f"{fname}.fa",
                'reference_seq': self._unitig_seqs[edge],
                'bed_content': bed_content

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


    def _custom_scoring_function(self, aligned_first, alignment_string, aligned_second,
                                first_to_ref, reference_to_first,
                                intersection_start, first_cl_dict, second_cl_dict,
                                commonSNPs, first_cl_start):
        info_printed = False
        """
        A custom distance scoring function for two sequences taking into account the artifacts of Flye consensus.
        alignment_string: a string consisting of '-', '.', '|' characters which correspond to indel, mismatch, match,
        respectively.
        Mismatches are penalized if they appear in commonSNPs.
        Indels are 1 point each if there are more than 5 of them in a contiguous block.
        Indel blocks start at the beginning or finish at the end are ignored.
        Variants that are covered by less than self._coverage_limits are assumed matches.
        """
        score = 0
        indel_length = 0
        indel_block_start = -1
        alignment_list = list(alignment_string)

        first_shift = intersection_start - first_cl_dict['start']
        second_shift = intersection_start - second_cl_dict['start']

        # read the contents of the bed.gz files that contain coverage information for each coordinate
        # list of lists with (start, end, coverage) for each coordinate interval
        cl1_bed_contents = first_cl_dict['bed_content']
        cl2_bed_contents = second_cl_dict['bed_content']


        for i, base in enumerate(alignment_list):
            if base not in "-.|":
                raise Exception("Unknown alignment sybmol!")

            # true coordinate = current coordinate on the alignment_string
            # + start of the intersection
            # - gaps in the target (or query) sequence thus far
            cl1_true_coor = first_shift + i - aligned_first[:i].count('-')
            cl2_true_coor = second_shift + i - aligned_second[:i].count('-')

            # ignore variants with low coverage
            if ((base == '-' or base == '.')
                    and (calculate_coverage(cl1_true_coor, cl1_bed_contents) < self._coverage_limit
                    or calculate_coverage(cl2_true_coor, cl2_bed_contents) < self._coverage_limit)):
                # TODO: this makes it a problem to calculate the true coordinate later on
                base = '|'

            if base == '-':
                # ignore the indels at the first or last position
                if i == 0:
                    indel_length = 1
                    indel_block_start = i
                    continue
                elif i == (len(alignment_list) - 1):
                    # trailing gaps are ignored
                    # break out of the loop without increasing the score
                    break
                else:
                    # start of an indel block if previous base was not a gap
                    if alignment_list[i - 1] != '-':
                        indel_block_start = i
                        indel_length = 1
                    # extend the indel block
                    else:
                        indel_length += 1

            # a contiguous gap ends if the previous base was a gap but not this one
            elif i != 0 and alignment_list[i - 1] == '-':
                if (indel_length >= self._indel_block_length_leniency and
                        indel_block_start != 0):
                    score += indel_length
                indel_length = 0

            # mismatch
            if base == '.':
                if info_printed == False:
                    info_printed = True
                    # logger.info(f"common SNP positions for the unitig: {commonSNPs}")

                mismatch_position = self._get_true_mismatch_position(
                    aligned_first,
                    first_to_ref,
                    reference_to_first,
                    i,
                    intersection_start - first_cl_start
                    ) + first_cl_start
                
                if mismatch_position in commonSNPs:
                    self._position_hit.value += 1
                    score += 1
                else:
                    self._position_miss.value += 1

        return score


    def _get_true_mismatch_position(self, cons_to_cons, cons_to_ref, reference, mismatch_index, first_cl_start):

        """ 
        The following three are outputs of alignments and contain gaps:
        cons_to_cons: consensus sequence aligned to another consensus sequence
        cons_to_ref: same consensus sequence aligned to reference
        reference: reference aligned to cons_to_ref
        mismatch_index: index of the mismatch in cons_to_cons
        """

        # TODO: this function may return a position that corresponds to a gap
        # in reference. Is this valid?

        # TODO: count mismatches too?
        # Find how many bases are there up to and including the mismatch_index
        true_pos_cons_to_cons = mismatch_index - cons_to_cons[:mismatch_index].count('-') + first_cl_start + 1
        
        # Find the index of base from cons_to_cons in cons_to_ref
        true_pos_cons_to_ref = -1
        bases = 0 # number of bases so far
        for i, b in enumerate(cons_to_ref):
            if b != '-':
                bases += 1
            if bases == true_pos_cons_to_cons:
                true_pos_cons_to_ref = i
                break

        # -1 indicates error
        if true_pos_cons_to_ref == -1:
            return true_pos_cons_to_ref
        
        # Find the true position of base matched with cons_to_ref in reference
        return true_pos_cons_to_ref - reference[:true_pos_cons_to_ref].count('-') + 1


    def cluster_distance_via_alignment(self, first_cl, second_cl, cl, edge, commonSNPs, debug=False):
        """
        Computes the distance between two clusters consensus'. The distance is based on the global alignment between the
        intersecting parts of the consensus'.
        first_cl: id (int)
        second_cl: id (int)
        cl: dataframe with columns 'read_name' and 'cluster' (id)
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

        reference_seq = first_cl_dict['reference_seq']
        aligned_first, aligned_second, edlib_aln = self._edlib_align(first_consensus_clipped, second_consensus_clipped)

        
        # check if alignment to reference is already computed
        cache_key = f"{edge}-{first_cl}-{first_cl_dict['start']}-{first_cl_dict['end']}"
        
        # Check if the result is already computed 
        already_computed = False
        with self._lock:
            if cache_key in self._alignment_cache:
                already_computed = True

        if already_computed:
            self._alignment_cache_hit.value += 1
            with self._lock:
                first_cl_to_ref, reference_aligned =  self._alignment_cache[cache_key]
        else:
            self._alignment_cache_miss.value += 1
            first_cl_to_ref, reference_aligned, _ = self._edlib_align(first_cl_dict['consensus'], reference_seq[first_cl_dict['start']:first_cl_dict['end']])
            # cache the reference alignment for re-use
            with self._lock:
                self._alignment_cache[cache_key] = [first_cl_to_ref, reference_aligned]

        edlib_score = self._custom_scoring_function(aligned_first, edlib_aln, aligned_second, first_cl_to_ref, reference_aligned, intersection_start,
                                                    first_cl_dict, second_cl_dict, commonSNPs, first_cl_dict['start'])
        
        # score is not normalized!
        return edlib_score




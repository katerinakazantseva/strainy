import pandas as pd
#pd.options.mode.chained_assignment = None
import pysam
import logging
from strainy.params import *

logger = logging.getLogger()


def build_adj_matrix(cl, data, SNP_pos, I, file, edge, R, only_with_common_snip=True):
    m = pd.DataFrame(-1, index=cl['ReadName'], columns=cl['ReadName'])
    logger.debug("Building adjacency matrix with " + str(m.shape[1]) + " reads")
    if only_with_common_snip==False:
        for i in range(1, m.shape[1]):
            first_read = m.index[i]
            for j in range(1, m.shape[1]):
                second_read=m.index[j]
                m[second_read][first_read] = distance(first_read, second_read, data, SNP_pos, R, only_with_common_snip=False)
    else:
        for i in range(1, m.shape[1]):
            first_read = m.index[i]
            bamfile = pysam.AlignmentFile(file, "rb")
            border1 = data[first_read]["Start"] + I
            border2 = data[first_read]["End"] - I
            if border2 <= 0:
                border2 = 1
            if border1 <= 0:
                border1 = 1

            for pos in [border1, border2]:
                for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools',
                                                   ignore_overlaps=False,
                                                   ignore_orphans=False,
                                                   truncate=True):
                    for pileupread in pileupcolumn.pileups:
                        second_read = pileupread.alignment.query_name

                        try:
                            if m[second_read][first_read] == -1:
                                m[second_read][first_read] = distance(first_read, second_read, data, SNP_pos, R)
                        except KeyError:
                            pass
    return m


def distance(read1, read2, data, SNP_pos, R, only_with_common_snip=True):
    d = -1
    firstSNPs = list(data[read1].keys())
    secondSNPs = list(data[read2].keys())
    keys=('Stop','Start')
    firstSNPs = [key for key in firstSNPs if key not in keys]
    secondSNPs= [key for key in secondSNPs if key not in keys]
    commonSNP = sorted(set(firstSNPs).intersection(secondSNPs).intersection(SNP_pos))
    #if len(commonSNP)>2 or (len(commonSNP) > 0 and only_with_common_snip == False) or len(SNP_pos)<=10: #???
    if 1==1:
        for snp in commonSNP:
            try:
                b1 = data[read1][snp]
                b2 = data[read2][snp]
                if b1 != b2 and len(b1) != 0 and len(b2) != 0:
                    if d == -1:
                        d = 0
                    d = d + 1
                elif b1 == b2:
                    if d == -1:
                        d = 0
                    d = d
            except:
                continue
            if d >= R:
                d = R
                break

            else:
                continue
    if len(commonSNP) == 0 and only_with_common_snip == False:
        intersect = set(range(data[read1]["Start"], data[read1]["End"])).intersection(
            set(range(data[read2]["Start"], data[read2]["End"])))
        if len(intersect) > 0:
            d = 0
        else:
            d = 1
    return d


def remove_edges(m, R):
    m_transformed = m
    m_transformed[m_transformed >= R] = -1
    return m_transformed


def change_w(m, R):
    m_transformed = m
    m_transformed[m_transformed == 0] = 0.001
    m_transformed[m_transformed == -1] = 0
    m_transformed[m_transformed >= R] = 0
    return m_transformed


def distance_clusters(edge,first_cl,second_cl, cons,cl, flye_consensus, only_with_common_snip=True):
    d = -1
    firstSNPs = list(cons[first_cl].keys())
    secondSNPs = list(cons[second_cl].keys())
    keys=('clSNP','clSNP2', 'Strange', 'Strange2','End','Start','Cov')
    firstSNPs = set([int(key) for key in firstSNPs if key not in keys])
    secondSNPs = set([int(key) for key in secondSNPs if key not in keys])
    commonSNP = set(sorted(firstSNPs.intersection(secondSNPs)))
    intersect = set(range(cons[first_cl]["Start"],cons[first_cl]["End"])).intersection(set(range(cons[second_cl]["Start"],cons[second_cl]["End"])))

    if only_with_common_snip == False and len(commonSNP) == 0 and len(intersect) > I:
        d = 0
    elif only_with_common_snip == True and len(set(cons[first_cl]["clSNP2"]).intersection(set(cons[second_cl]["clSNP2"]))) == 0:
        d = 1
    elif len(intersect) > I:
        d = flye_consensus.cluster_distance_via_alignment(first_cl, second_cl, cl, edge, commonSNP)
    else:
        d = 1
    return d

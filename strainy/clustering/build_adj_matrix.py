import logging

import pandas as pd
from scipy.spatial.distance import cdist

from strainy.params import *

logger = logging.getLogger()
#pd.options.mode.chained_assignment = None

class DistanceWrapper():
    # Wrapper for calling cdist with custom distance function
    def __init__(self, cl, data, SNP_pos, R, only_with_common_snip):
        self.cl = cl
        self.data = data
        self.SNP_pos = SNP_pos
        self.R = R
        self.only_with_common_snip = only_with_common_snip

    def distance_wrapper(self, first_read, second_read):
        return distance(first_read[0],
                        second_read[0],
                        self.data,
                        self.SNP_pos,
                        self.R,
                        self.only_with_common_snip)


def build_adj_matrix(cl, data, SNP_pos, I, file, edge, R, only_with_common_snip=True):
    m = pd.DataFrame(-1.0, index=cl['ReadName'], columns=cl['ReadName'])
    logger.debug("Building adjacency matrix with " + str(m.shape[1]) + " reads")
    if only_with_common_snip==False:
        dw = DistanceWrapper(cl, data, SNP_pos, R, only_with_common_snip)
        result = cdist(cl['ReadName'].to_frame(), cl['ReadName'].to_frame(), dw.distance_wrapper)

        # Set the first row and the column to -1
        try:
            result[0,:] = -1
            result[:,0] = -1
        except IndexError:
            pass

        result_df = pd.DataFrame(result, 
                         index=cl['ReadName'],
                         columns=cl['ReadName'])
    else:

        dw = DistanceWrapper(cl, data, SNP_pos, R, only_with_common_snip)
        result = cdist(cl['ReadName'].to_frame(), cl['ReadName'].to_frame(), dw.distance_wrapper)

        result[0,:] = -1
        result_df = pd.DataFrame(result, 
                        index=cl['ReadName'],
                        columns=cl['ReadName'])

    return result_df


def distance(read1, read2, data, SNP_pos, R, only_with_common_snip=True):
    d = -1
    firstSNPs = list(data[read1].keys())
    secondSNPs = list(data[read2].keys())
    keys=('End','Start', 'Rclip', 'Lclip')
    firstSNPs = [key for key in firstSNPs if key not in keys]
    secondSNPs= [key for key in secondSNPs if key not in keys]
    commonSNP = sorted(set(firstSNPs).intersection(secondSNPs).intersection(SNP_pos))

    if read1 == read2:
        return 0

    if only_with_common_snip:
        intersect = max(min(data[read1]["End"], data[read2]["End"]) - max(data[read1]["Start"], data[read2]["Start"]), 0)
        if intersect < I:
            return -1.0

        if len(commonSNP) == 0:
            return -1.0

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
            except:
                continue

        d = d / intersect

    if len(commonSNP) == 0 and only_with_common_snip == False:
        intersect = max(min(data[read1]["End"], data[read2]["End"]) - max(data[read1]["Start"], data[read2]["Start"]), 0)
        if intersect > 0:
            d = 0
        else:
            d = 1

    return float(d)


def remove_edges(m, R):
    m_transformed = m
    m_transformed[m_transformed > R] = -1
    return m_transformed


def change_w(m, R):
    m_transformed = m
    m_transformed[m_transformed == 0] = -10
    m_transformed[m_transformed == -1] = 0
    m_transformed[m_transformed > R] = 0
    m_transformed[m_transformed == -10] = 0.000001

    return m_transformed


def distance_clusters(edge,first_cl,second_cl, cons,cl, flye_consensus, only_with_common_snip=True):
    d = -1
    firstSNPs = list(cons[first_cl].keys())
    secondSNPs = list(cons[second_cl].keys())
    keys=('clSNP','clSNP2', 'Strange', 'Strange2','End','Start','Cov')
    firstSNPs = set([int(key) for key in firstSNPs if key not in keys])
    secondSNPs = set([int(key) for key in secondSNPs if key not in keys])
    commonSNP = set(sorted(firstSNPs.intersection(secondSNPs)))
    intersect = max(min(cons[first_cl]["End"], cons[second_cl]["End"]) - max(cons[first_cl]["Start"], cons[second_cl]["Start"]), 0)

    if only_with_common_snip == False and len(commonSNP) == 0 and intersect > I:
        d = 0
    elif only_with_common_snip == True and len(set(cons[first_cl]["clSNP2"]).intersection(set(cons[second_cl]["clSNP2"]))) == 0:
    #elif only_with_common_snip == True and len(commonSNP) == 0:
        d = 1
    elif intersect > I:
        d = (flye_consensus.cluster_distance_via_alignment(first_cl, second_cl, cl, edge, commonSNP))/intersect
    else:
        d = 1
    return float(d)

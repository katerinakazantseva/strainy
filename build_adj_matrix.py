import pandas as pd
import pysam


def build_adj_matrix2(cl, data, SNP_pos, I, file, edge, R, only_with_common_snip=True):
    m = pd.DataFrame(-1, index=cl['ReadName'], columns=cl['ReadName'])
    if only_with_common_snip==False:
        for i in range(1, m.shape[1]):
            first_read = m.index[i]
            for j in range(1, m.shape[1]):
                second_read=m.index[j]
                m[second_read][first_read] = distance(first_read, second_read, data, SNP_pos, R, only_with_common_snip=False)

    else:
        for i in range(1, m.shape[1]):
            print(str(i) + "/" + str(m.shape[1]) + " Reads processed \r", end="")
            first_read = m.index[i]
            bamfile = pysam.AlignmentFile(file, "rb")
            border1 = data[first_read]["Start"] + I
            border2 = data[first_read]["Stop"] - I
            if border2 <= 0:
                border2 = 1
            if border1 <= 0:
                border1 = 1

            for pos in [border1, border2]:
                for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools',
                                               ignore_overlaps=False,
                                               ignore_orphans=False, truncate=True):
                    for pileupread in pileupcolumn.pileups:
                        second_read = pileupread.alignment.query_name

                        try:
                            if m[second_read][first_read] == -1:
                                m[second_read][first_read] = distance(first_read, second_read, data, SNP_pos, R)
                        except:
                            KeyError
    return (m)


def distance(read1, read2, data, SNP_pos, R, only_with_common_snip=True):
    d = -1
    firstSNPs = list(data[read1].keys())
    firstSNPs.remove('Start')
    firstSNPs.remove('Stop')
    secondSNPs = list(data[read2].keys())
    secondSNPs.remove('Start')
    secondSNPs.remove('Stop')
    commonSNP = sorted(set(firstSNPs).intersection(secondSNPs).intersection(SNP_pos))
    # for snp in SNP_pos:
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
        intersect = set(range(data[read1]["Start"], data[read1]["Stop"])).intersection(
            set(range(data[read2]["Start"], data[read2]["Stop"])))
        if len(intersect) > 0:
            d = 0
        else:
            d = 1
    return (d)


def remove_edges(m, R):
    m_transformed = m
    m_transformed[m_transformed >= R] = -1
    return (m_transformed)


def change_w(m, R):
    m_transformed = m
    m_transformed[m_transformed == 0] = 0.001
    m_transformed[m_transformed == -1] = 0
    m_transformed[m_transformed >= R] = 0
    return (m_transformed)

import pandas as pd
import pysam


def build_adj_matrix2 (cl,data,SNP_pos,I, file, edge, R):

    m = pd.DataFrame(-1, index=cl['ReadName'], columns=cl['ReadName'])

    for i in range(1,m.shape[1]):
        print(str(i)+"/"+str(m.shape[1])+" Reads processed \r", end="")
        first_read=m.index[i]
        bamfile = pysam.AlignmentFile(file, "rb")
        border1=data[first_read]["Start"] + I
        border2=data[first_read]["Stop"] - I
        if border2<=0:
            border2=1
        if border1<=0:
            border1=1

        for pos in [border1, border2]:
            for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools',
                                               ignore_overlaps=False,
                                               ignore_orphans=False, truncate=True):
                for pileupread in pileupcolumn.pileups:
                    second_read = pileupread.alignment.query_name

                    try:
                        if m[second_read][first_read]==-1:
                            m[second_read][first_read] = distance(first_read, second_read, data, SNP_pos, R)
                    except: KeyError
    return (m)


def distance(read1,read2,data, SNP_pos, R):
    d=-1
    for snp in SNP_pos:
        try:
            b1=data[read1][snp]
            b2=data[read2][snp]
            if b1 != b2 and len(b1)!=0 and  len(b2)!=0:
                if d==-1:
                    d=0
                d=d+1
            elif b1 == b2:
                if d==-1:
                    d=0
                d=d
        except:
            continue
        if d>=R:
            d=R
            break

        else:
            continue
    return (d)

def remove_edges (m, R):
    m_transformed=m
    m_transformed[m_transformed >= R] = -1
    return (m_transformed)

def change_w (m, R):
    m_transformed = m
    m_transformed[m_transformed ==0] = 0.001
    m_transformed[m_transformed == -1] = 0
    m_transformed[m_transformed >=R] = 0
    return (m_transformed)

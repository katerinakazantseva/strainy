import subprocess
import pysam
import os
import re
from collections import Counter

os.environ["PATH"] += os.pathsep + "/usr/local/bin"


def read_snp(snp,edge, bam, AF):
    SNP_pos = []
    if snp==None:
        snpos = 'bcftools mpileup -r {} {} --no-reference -I --no-version --annotate FORMAT/AD | bcftools query -f  "%CHROM %POS [ %AD %DP]\n" >output/vcf/vcf_{}.txt'.format(edge,bam,edge)
        subprocess.check_output(snpos, shell=True, capture_output=False)
        with open("output/vcf/vcf_%s.txt" % edge) as f:
            lines = f.readlines()
            for line in lines:
                try:
                    AlFreq = int(str(line.split()[2]).split(',')[2])/int(line.split()[3])
                except IndexError:
                    AlFreq = 0
                if AlFreq>AF:
                    SNP_pos.append(line.split()[1])
    else:
        vcf = open(snp, "rt")
        for line in vcf:
            if line.split()[0] == edge:
                SNP_pos.append(line.split()[1])
    print(str(len(SNP_pos)) + " SNPs found")
    return SNP_pos


def read_bam(bam, edge, SNP_pos, clipp, min_mapping_quality, min_al_len, de_max):
    bamfile = pysam.AlignmentFile(bam, "rb")
    data = {}
    ln = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4]
    for read in bamfile.fetch(edge):

        clipping = 0
        start = read.get_reference_positions()[0]
        stop = read.get_reference_positions()[len(read.get_reference_positions()) - 1]
        de=float(re.sub(".*de',\s", "", str(str(read).split('\t')[11]), count=0, flags=0).split(')')[0])
        for i in read.cigartuples:
            if i[0] == 4 or i[0] == 5:
                if i[1] > clipp:
                    clipping = 1

        if read.mapping_quality>min_mapping_quality and de < de_max and (((clipping == 0 and (stop - start) > min_al_len) and (
                int(start) != 0 and int(stop) != int(ln) - 1)) or int(start) == 0 or int(stop) == int(ln) - 1):
            data[read.query_name] = {}
            data[read.query_name]["Start"]=start
            data[read.query_name]["Stop"] = stop
            tags = read.get_tags()[9]

            left_clip = []
            right_clip = []

            if re.search("SA", str(tags)):
                if start == 0:
                    for i in tags[1].split(';'):
                        if len(i)>0:
                         left_clip.append(i.split(',')[0])
                if stop ==  int(ln) - 1:
                    for i in tags[1].split(';'):
                        if len(i) > 0:
                            right_clip.append(i.split(',')[0])
            data[read.query_name]["Rclip"]=list(set(right_clip))
            data[read.query_name]["Lclip"]=list(set(left_clip))
        else:
            continue

    for pos in SNP_pos:
        for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools', min_base_quality=0,
                                           ignore_overlaps=False, min_mapping_quality=min_mapping_quality,
                                           ignore_orphans=False, truncate=True):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    try:
                        data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[pileupread.query_position]
                    except KeyError:
                        continue
    bamfile.close()
    return data


def build_data_cons(cl,SNP_pos, data):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
    cons = {}
    for cluster in clusters:
        cons=cluster_consensuns(cl,cluster,SNP_pos, data, cons)
    return cons


def cluster_consensuns(cl, cluster, SNP_pos, data, cons):
    strange = 0
    val = {}
    clSNP = []
    for pos in SNP_pos:
        npos = []
        for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
            try:
                npos.append(data[read][pos])
            except KeyError:
                continue
        try:
            if len(npos) >= 2:
                if int(Counter(npos).most_common()[0][1]) >= 2:
                    val[pos] = Counter(npos).most_common()[0][0]
            if int(Counter(npos).most_common()[1][1]) >= 2:
                strange = 1
                clSNP.append(pos)
        except(IndexError):
            continue

    if strange == 1:
        val["Strange"] = 1
    else:
        val["Strange"] = 0
    val["clSNP"] = clSNP

    clStart = 1000000000000  # change fo ln
    clStop = 0
    clCov=0

    for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
        start = int(data[read]["Start"])
        stop = int(data[read]["Stop"])
        clCov = clCov+(stop-start)
        if start < clStart:
            clStart = start
        if stop > clStop:
            clStop = stop

    clCov=clCov/(clStop-clStart)
    val["Stop"] = clStop
    val["Start"] = clStart
    val["Cov"]=clCov
    cons[cluster] = val
    # val = {
    # "Stop": end of the last read,
    # "Start: beginning of the first read,
    # "Cov: X cluster coverage,
    # "clSNP = [], positions where strange happens
    # pos: most common,
    # pos: most common,
    # ...
    # }
    return (cons)

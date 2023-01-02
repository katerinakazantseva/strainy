import subprocess
import pysam
import os
import re
from collections import Counter

from metaphase.params import *

os.environ["PATH"] += os.pathsep + "/usr/local/bin"


def read_snp(snp, edge, bam, AF,cluster=None):
    SNP_pos = []
    if snp==None:
        if cluster==None:
            snpos = ('bcftools mpileup -r {} {} --no-reference -I --no-version --annotate FORMAT/AD 2>/dev/null ' +
                     '| bcftools query -f  "%CHROM %POS [ %AD %DP]\n" >{}/vcf/vcf_{}.txt').format(edge, bam, MetaPhaseArgs.output, edge)
            subprocess.check_output(snpos, shell=True, capture_output=False)
            #subprocess.call(snpos, shell=True, stderr=subprocess.DEVNULL)
            with open("%s/vcf/vcf_%s.txt" % (MetaPhaseArgs.output, edge)) as f:
                lines = f.readlines()
                for line in lines:
                    try:
                        AlFreq = int(str(line.split()[2]).split(',')[2]) / int(line.split()[3])
                    except(IndexError):
                        AlFreq = 0
                    if AlFreq > AF:
                        SNP_pos.append(line.split()[1])
        else:
            snpos = ('bcftools mpileup -f {} -r {} {}  -I --no-version --annotate FORMAT/AD 2>/dev/null ' +
                     '| bcftools query -f  "%CHROM %POS %ALT [ %AD %DP]\n" >{}/vcf/vcf_{}_{}.txt').format(MetaPhaseArgs.fa,edge, bam, MetaPhaseArgs.output, edge, cluster)
            subprocess.check_output(snpos, shell=True, capture_output=False)
            with open("%s/vcf/vcf_%s_%s.txt" % (MetaPhaseArgs.output, edge, cluster)) as f:
                lines = f.readlines()
                for line in lines:
                    try:
                        AlFreq = int(str(line.split()[3]).split(',')[2]) / int(line.split()[4])
                    except(IndexError):
                        AlFreq = 0
                    if AlFreq > AF:
                        SNP_pos.append(line.split()[1])
                    if int(str(line.split()[3]).split(',')[0])==0 and int(str(line.split()[3]).split(',')[1])>min_reads_cluster and len(str(line.split()[2]).split(','))>1:
                        SNP_pos.append(line.split()[1])

    else:
        vcf = open(snp, "rt")
        for line in vcf:
            if line.split()[0] == edge:
                SNP_pos.append(line.split()[1])
    # print(str(len(SNP_pos)) + " SNPs found")
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
        if read.mapping_quality>=min_mapping_quality  and de < de_max and (((clipping == 0 and (stop - start) > min_al_len) and (
                int(start) != 0 and int(stop) != int(ln) - 1)) or int(start) < 5  or int(stop) > int(ln) - 5):

            data[read.query_name] = {}
            data[read.query_name]["Start"]=start
            data[read.query_name]["Stop"] = stop
            try:
                tags = read.get_tags()[9]
            except (IndexError):
                tags=""
            left_clip = {}
            right_clip = {}

            if re.search("SA", str(tags)):
                if read.cigartuples[0][0]==4 and read.cigartuples[len(read.cigartuples)-1][0]!=4:
                    for i in tags[1].split(';'):
                        if len(i)>0 and int(i.split(',')[4])>20:
                         orient = '+' if read.is_reverse==False else '-'
                         left_clip[i.split(',')[0]]=[i.split(',')[2], orient]
                         break
                if read.cigartuples[0][0] != 4 and read.cigartuples[len(read.cigartuples) - 1][0] == 4:
                    for i in tags[1].split(';'):
                        if len(i) > 0 and int(i.split(',')[4])>20:
                            orient = '+' if read.is_reverse == False else '-'
                            right_clip[i.split(',')[0]] = [i.split(',')[2], orient]
                else:
                    for i in tags[1].split(';'):
                        if len(i) > 0 and int(i.split(',')[4]) > 20:
                            l=int(re.match('.*?([0-9]+)$', re.sub(r"M.*$", "", str(i.split(',')[3]))).group(1))
                            if l==read.cigartuples[0][1]:
                                orient = '+' if read.is_reverse == False else '-'
                                left_clip[i.split(',')[0]] = [i.split(',')[2], orient]
                            elif l==read.cigartuples[len(read.cigartuples) - 1][1]:
                                orient = '+' if read.is_reverse == False else '-'
                                right_clip[i.split(',')[0]] = [i.split(',')[2], orient]
                            elif l>read.cigartuples[0][1]:
                                orient = '+' if read.is_reverse == False else '-'
                                right_clip[i.split(',')[0]] = [i.split(',')[2], orient]
                            elif l> read.cigartuples[len(read.cigartuples) - 1][1]:
                                orient = '+' if read.is_reverse == False else '-'
                                left_clip[i.split(',')[0]] = [i.split(',')[2], orient]
            data[read.query_name]["Rclip"]=right_clip
            data[read.query_name]["Lclip"]=left_clip
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

                    except (KeyError):
                        continue
    bamfile.close()
    return (data)



def build_data_cons(cl, SNP_pos, data,edge):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
    cons = {}
    for cluster in clusters:
        cons = cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge)
    return cons


def cluster_consensuns(cl, cluster, SNP_pos, data, cons, edge):
    strange = 0
    strange2 = 0
    val = {}
    clSNP = []
    clSNP2 = []

    reads = list(cl.loc[cl['Cluster'] == cluster, 'ReadName'])
    with open('%s/clusters/reads_%s_%s.txt' % (MetaPhaseArgs.output, edge, cluster), 'w') as fp:
        for line in reads:
            fp.write(str(line))
            fp.write("\n")

    #Create bam for cluster
    pysam.samtools.view("-N", '%s/clusters/reads_%s_%s.txt' % (MetaPhaseArgs.output, edge, cluster),
                        "-o", '%s/bam/clusters/%s_%s.bam' % (MetaPhaseArgs.output, edge, cluster), MetaPhaseArgs.bam, edge,catch_stdout=False)
    pysam.samtools.index('%s/bam/clusters/%s_%s.bam' % (MetaPhaseArgs.output, edge, cluster))

    clSNP2 = read_snp(MetaPhaseArgs.snp, edge, '%s/bam/clusters/%s_%s.bam' % (MetaPhaseArgs.output, edge, cluster), AF, cluster)

    try:
        if len(clSNP2) > 0 and max([int(clSNP2[i + 1]) - int(clSNP2[i])
                                    for i in range(0, len(clSNP2) - 1)]) > 1.5 * 1000:
            strange2 = 1

    except(ValueError):
        pass

    for pos in SNP_pos:
        npos = []
        for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
            try:
                npos.append(data[read][pos])
            except(KeyError):
                continue

        try:
            if len(npos) > 2:
                if int(Counter(npos).most_common()[0][1]) > 2:
                    val[pos] = Counter(npos).most_common()[0][0]
            if int(Counter(npos).most_common()[1][1]) >= unseparated_cluster_min_reads:
                strange = 1
                clSNP.append(pos)

        except(IndexError):
            continue

    val["clSNP"] = clSNP
    val["clSNP2"] = clSNP2

    clStart = 1000000000000  # change fo ln
    clStop = 0
    clCov = 0

    for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
        try:
            start = int(data[read]["Start"])
            stop = int(data[read]["Stop"])
            clCov = clCov + (stop - start)
            if start < clStart:
                clStart = start
            if stop > clStop:
                clStop = stop

        except(KeyError):
            pass

    try:
        if (int(clSNP2[0]) - int(clStart)) > 1.5 * I or \
                int(clStop) - int(clSNP2[len(clSNP2) - 1]) > 1.5 * I:
            strange2 = 1
    except(IndexError):
        pass

    val["Strange"] = int(strange == 1)
    val["Strange2"] = int(strange2 == 1)
    clCov = clCov / (clStop - clStart)
    val["Stop"] = clStop
    val["Start"] = clStart
    val["Cov"] = clCov
    cons[cluster] = val

    return cons

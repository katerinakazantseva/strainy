import subprocess
import pysam
from Bio import SeqIO
import re
from collections import Counter


def read_snp(snp,edge, bam, AF):
    SNP_pos = []
    if snp==None:
        snpos = 'bcftools mpileup -r {} {} --no-reference -I --no-version --annotate FORMAT/AD | bcftools query -f  "%CHROM %POS [ %AD %DP]\n" >output/vcf/vcf_{}.txt'.format(edge,bam,edge)
        subprocess.check_output(snpos, shell=True, capture_output=False)
        with open("output/vcf/vcf_%s.txt" % edge) as f:
            lines = f.readlines()
            for line in lines:
                try:
                    AlFreq=int(str(line.split()[2]).split(',')[2])/int(line.split()[3])
                except(IndexError):  AlFreq=0
                if AlFreq>AF:
                    SNP_pos.append(line.split()[1])

    else:
        vcf = open(snp, "rt")
        for line in vcf:
            if line.split()[0] == edge:
                SNP_pos.append(line.split()[1])
    print(str(len(SNP_pos)) + " SNPs found")
    return(SNP_pos)


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
        #and read.is_supplementary == False
        if read.mapping_quality>=min_mapping_quality  and de < de_max and (((clipping == 0 and (stop - start) > min_al_len) and (
                int(start) != 0 and int(stop) != int(ln) - 1)) or int(start) < 3  or int(stop) == int(ln) - 1):

            data[read.query_name] = {}
            #print(read.query_name)
            data[read.query_name]["Start"]=start
            data[read.query_name]["Stop"] = stop
            try:
                tags = read.get_tags()[9]
            except (IndexError):
                tags=""

            left_clip = {}
            right_clip = {}

            if re.search("SA", str(tags)):

                #if start==0 and stop ==  int(ln) - 1:
                #if read.cigartuples[0][0] == 4 and read.cigartuples[len(read.cigartuples) - 1][0] == 4:
                    #print(read.query_name)

                    # print(read.query_alignment_end-read.query_alignment_length)
                    # print(read.query_alignment_end)
                    #print(tags)
                    #print(read.cigartuples)
                    #print()
                #if start == 0:
                #print(read.query_name)
                #fields = read.to_string().split()
                #read_id, flags, chr_id, position = fields[0:4]
                #ref_id, ref_start = fields[2], int(fields[3])
                #print(ref_id)
                #print(ref_start)
                if read.cigartuples[0][0]==4 and read.cigartuples[len(read.cigartuples)-1][0]!=4:
                    for i in tags[1].split(';'):
                        if len(i)>0 and int(i.split(',')[4])>20:
                         #print(i)
                         orient = '+' if read.is_reverse==False else '-'
                         left_clip[i.split(',')[0]]=[i.split(',')[2], orient]
                         break
                #if stop ==  int(ln) - 1:
                if read.cigartuples[0][0] != 4 and read.cigartuples[len(read.cigartuples) - 1][0] == 4:
                    for i in tags[1].split(';'):
                        if len(i) > 0 and int(i.split(',')[4])>20:
                            #print(i)
                            orient = '+' if read.is_reverse == False else '-'
                            right_clip[i.split(',')[0]] = [i.split(',')[2], orient]
                else:
                    #print(read.cigartuples[0][1])
                    #print(read.cigartuples[len(read.cigartuples) - 1][1])
                    for i in tags[1].split(';'):
                        if len(i) > 0 and int(i.split(',')[4]) > 20:
                            l=int(re.match('.*?([0-9]+)$', re.sub(r"M.*$", "", str(i.split(',')[3]))).group(1))
                            if l==read.cigartuples[0][1]:
                                orient = '+' if read.is_reverse == False else '-'
                                left_clip[i.split(',')[0]] = [i.split(',')[2], orient]
                            elif l==read.cigartuples[len(read.cigartuples) - 1][1]:
                                orient = '+' if read.is_reverse == False else '-'
                                right_clip[i.split(',')[0]] = [i.split(',')[2], orient]



                #print()
            #print(left_clip)


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




def build_data_cons(cl,SNP_pos, data):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
    cons = {}
    for cluster in clusters:
        cons=cluster_consensuns(cl,cluster,SNP_pos, data, cons)
    return(cons)


def cluster_consensuns(cl,cluster,SNP_pos, data, cons):
    strange = 0
    val = {}
    clSNP = []
    for pos in SNP_pos:
        npos = []
        for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
            try:
                npos.append(data[read][pos])
            except(KeyError):
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
        try:
            start = int(data[read]["Start"])
            stop = int(data[read]["Stop"])
            clCov=clCov+(stop-start)
            if start < clStart:
                clStart = start
            if stop > clStop:
                clStop = stop
        except(KeyError):
            pass

    clCov=clCov/(clStop-clStart)
    val["Stop"] = clStop
    val["Start"] = clStart
    val["Cov"]=clCov
    cons[cluster] = val
    return (cons)

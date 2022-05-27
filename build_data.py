import subprocess
import pysam
from Bio import SeqIO
import re


def read_snp(snp,edge, bam, AF):
    SNP_pos = []
    if snp==None:
        snpos = 'bcftools mpileup -r {} {} --no-reference -I --no-version --annotate FORMAT/AD | bcftools query -f  "%CHROM %POS [ %AD %DP]\n" >output/vcf/vcf_{}.txt'.format(edge,bam,edge)
        #snpos = 'bcftools mpileup -r {}  /Users/ekaterina.kazantseva/MT/test_data/SRR13128014_labeled.bam   --no-reference -I --no-version --annotate FORMAT/AD | bcftools query -f  "%CHROM %POS [ %AD %DP]\n" >output/vcf/vcf_{}.txt'.format(edge,edge)
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
    #print(SNP_pos)
    return(SNP_pos)




def read_bam_old(bam,edge,SNP_pos,clipp,min_mapping_quality,min_al_len,de_max):
    bamfile = pysam.AlignmentFile(bam, "rb")
    data = {}
    ln = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4]
    for pos in SNP_pos:
    #for pos in range(1,int(ln)+1, 500):
        #print(pos)
        for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools', min_base_quality=0,
                                           ignore_overlaps=False,min_mapping_quality=min_mapping_quality,
                                           ignore_orphans=False, truncate=True):
            for pileupread in pileupcolumn.pileups:

                clipping=0
                start = pileupread.alignment.get_reference_positions()[0]
                stop = pileupread.alignment.get_reference_positions()[
                    len(pileupread.alignment.get_reference_positions()) - 1]

                de=float(str(pileupread).split('\t')[11].split('), (')[8].split(',')[1])
                for i in pileupread.alignment.cigartuples:
                    if i[0] == 4 or i[0] == 5:
                        if i[1]>clipp:
                            clipping=1



                if (clipping==1 or (stop-start)<min_al_len or de>de_max) and (int(start)!=0 and int(stop)!=int(ln)-1):
                    continue



                else:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        try:
                            data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[
                                pileupread.query_position]

                        except (KeyError):
                            data[pileupread.alignment.query_name] = {}
                            data[pileupread.alignment.query_name]["Start"] = pileupread.alignment.get_reference_positions()[
                            0]
                            data[pileupread.alignment.query_name]["Stop"] = pileupread.alignment.get_reference_positions()[
                                len(pileupread.alignment.get_reference_positions()) - 1]

                            data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[
                             pileupread.query_position]

    bamfile.close()
    return(data)


def read_bam(bam, edge, SNP_pos, clipp, min_mapping_quality, min_al_len, de_max):
    bamfile = pysam.AlignmentFile(bam, "rb")
    data = {}
    ln = pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4]

    for read in bamfile.fetch(edge):
        clipping = 0
        start = read.get_reference_positions()[0]
        stop = read.get_reference_positions()[len(read.get_reference_positions()) - 1]
        #de = float(str(read).split('\t')[11].split('), (')[7].split(',')[1])
        de=float(re.sub(".*de',\s", "", str(str(read).split('\t')[11]), count=0, flags=0).split(')')[0])
        #print(de)
        for i in read.cigartuples:
            if i[0] == 4 or i[0] == 5:
                if i[1] > clipp:
                    clipping = 1
        if read.mapping_quality>min_mapping_quality and de < de_max and (((clipping == 0 and (stop - start) > min_al_len) and (
                int(start) != 0 and int(stop) != int(ln) - 1)) or int(start) == 0  or int(stop) == int(ln) - 1):
            data[read.query_name] = {}
            data[read.query_name]["Start"]=start
            data[read.query_name]["Stop"] = stop
        else:
            #print(read.query_name)
            #print(read.mapping_quality)
            #print(start)
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
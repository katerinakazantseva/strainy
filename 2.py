import csv
import pysam
from Bio import SeqIO
import pandas as pd
import numpy as np

bam="/Users/ekaterina.kazantseva/MT/test_data/test.bam"
snp="/Users/ekaterina.kazantseva/MT/test2.vcf"
edge='edge_8'
R=1
#READS SNP
print ("### Reading SNPs...")
SNP_pos=[]
vcf=open(snp, "rt")
for line in vcf:
    if line.split()[0]==edge:
      SNP_pos.append(line.split()[1])

#READ READS AND POSITIONS
print ("### Reading Reads...")

RName=[]
Seq=[]
Pos=[]

bamfile = pysam.AlignmentFile(bam, "rb")
iter = bamfile.fetch(edge)



for pos in SNP_pos:
    #print("new snip")
    for pileupcolumn in bamfile.pileup("edge_8", int(pos)-1, int(pos), stepper='samtools', min_base_quality=0, ignore_overlaps=False,
                                   ignore_orphans=False, truncate=True):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:

                RName.append(pileupread.alignment.query_name)
                Seq.append(pileupread.alignment.query_sequence[pileupread.query_position])
                Pos.append(pos)
            #print('\tbase in read %s = %s' %
                 # (pileupread.alignment.query_name,
                  # pileupread.alignment.query_sequence[pileupread.query_position]))

bamfile.close()

d = {'ReadName': RName, 'Seq': Seq, 'Pos': Pos}
df = pd.DataFrame(data=d)

df = df.drop(df[(df['Seq'] == 'None')].index).reset_index()
cl=pd.DataFrame(data={'ReadName': RName})
cl=cl.drop_duplicates(subset=['ReadName'])
cl['Cluster'] = 'NA'


print(cl)

#for x in iter:

 #   if int(str(x).split()[3])!=0:
  #      N1=0
   #     N2=0
    #    print(str(x).split()[0])

     #   for i in x.cigartuples:
      #      if i[0]==1:
       #      N1=N1+i[1]
        #    if i[0]==2:
         #    N2=N2+i[1]
        ##print(x.reference_start)
        #print(N1)
        #print(N2)
        #print(x.reference_start-N1+N2)
        #Pos.append(x.reference_start-N1+N2)
        #Pos.append(x.reference_start - N1 + N2+1)
        #try:
           # print(str(x).split()[9][8108-(x.reference_start - N1 + N2+1):8108-(x.reference_start - N1 + N2+1)+10])
           # print(x.reference_start)
           # print(x.cigartuples)
           #print(N1)
           # print(N2)
        #except:
           # continue
        #RName.append(str(x).split()[0])
        #Seq.append(str(x).split()[9])
    #else:
     #   continue
        #print(str(x).split()[0])
        #print(x.query_alignment_start)
        #print(x.cigartuples)
        #N1=0
        #N2=0
        #for i in x.cigartuples:
         #   if i[0]==1:
          #   N1=N1+i[1]
           # if i[0]==2:
            # N2=N2+i[1]
        #print(N1)
        #print(N2)
        #Pos.append (-int(x.query_alignment_start))




#for read in bamfile.fetch(edge, start=8108, end=8109):
    #try:
    #alignedRefPositions = read.get_reference_positions()
    #refStart = alignedRefPositions[0]
    #refSequence = read.get_reference_sequence()
    #readSequence = read.query_alignment_sequence
     #   print(read.query_name)

    #infer_query_length
      #  print(read.query_alignment_start)
        #print(read.infer_query_length()) длина рида
        #print(read.cigartuples)
       # read("  ")
    #except:
     #   continue


#print(df['Pos'])





#CALCULATE DISTANCE and BUILD graph
print ("### Calculatind distances/Building adj matrix...")


def old_distance(read1,read2,pos1,pos2):
    d=-1
    for i in SNP_pos:
        #print(i)
        if int(i)>=int(pos1) and int(i)>=int(pos2):

            try:
                print(read1[int(i)-int(pos1):int(i)-int(pos1)+5])
                print(read2[int(i)-int(pos2):int(i)-int(pos2)+5])
                if read1[int(i)-int(pos1)]!=read2[int(i)-int(pos2)]:
                    if d==-1:
                        d=0
                    d=d+1
                else:
                    d=0
            except (IndexError):
                continue
        if d>=R:
            d=R
            break

        else:
            continue



    return (d)

print (df)
read=df.loc[ 0,'ReadName']
#print (df[]['Seq'].values)


#print(df.loc[df['ReadName'] == read][df['Pos'] != '8108']['Seq'].values)

#

def distance(read1,read2,df):
    d=-1
    i=0
    for snp in SNP_pos:
        i=i+1
        try:
            b1=df.loc[( df.ReadName == read1) & (df.Pos == '%s' %snp)]['Seq'].values
            b2 = df.loc[(df.ReadName == read2) & (df.Pos == '%s' % snp)]['Seq'].values
            #print(len(b1))
            if  len(b1)==0 or  len(b2)==0 or b1 != b2:
                #print("not same")
                if d==-1:
                    d=0
                d=d+1
            else:
                if d==-1:
                    d=0
                #print("same")

        except:
            print("out")
            continue
        if d>=R:
            d=R
            break

        else:
            continue



    return (d)
#read1=df.loc[ 112,'ReadName']
#read2=df.loc[ 113,'ReadName']
#pos1=df.loc[ 112,'Pos']+1
#pos2=df.loc[ 113,'Pos']+2

#print(df.loc[ 112,'ReadName'])

#print(distance(read1, read2,df))



def old_build_adj_matrix (df):
    m = pd.DataFrame(-1, index=df['ReadName'], columns=df['ReadName'])
    #m = pd.DataFrame(0, index=[1,2,3,4],columns=[1,2,3,4])
    for i in range(1,m.shape[1]):
        print(str(i)+"/"+str(m.shape[1])+" Reads processed \r" , end="")
        first_read=m.index[i]
        for j in range(0,i):
            second_read = m.index[j]
            m[second_read][first_read] = distance(str(df[(df['ReadName'] ==first_read)]['Seq'].values).strip('[]').strip("'"), str(df[(df['ReadName'] ==second_read)]['Seq'].values).strip('[]').strip("'"),str(df[(df['ReadName'] ==first_read)]['Pos'].values).strip('[]').strip("'"), str(df[(df['ReadName'] ==second_read)]['Pos'].values).strip('[]').strip("'"))

    return (m)


def build_adj_matrix (cl,df):
    m = pd.DataFrame(-1, index=cl['ReadName'], columns=cl['ReadName'])
    #m = pd.DataFrame(0, index=[1,2,3,4],columns=[1,2,3,4])
    for i in range(1,m.shape[1]):
        print(str(i)+"/"+str(m.shape[1])+" Reads processed \r" , end="")
        first_read=m.index[i]
        for j in range(0,i):
            second_read = m.index[j]
            m[second_read][first_read] = distance(first_read,second_read, df)

    return (m)


m=build_adj_matrix (cl,df)
#m.to_csv("adj_M_%s.csv" % edge)

#m=pd.read_csv("adj_M_edge_8.csv")
print(m)


print("### Removing overweighed egdes...")

def remove_edges (m, R):
    m[m >= R] = -1
    return (m)

remove_edges (m, R)

print(m)

print("### Searching connected components...")


#m=m.set_index('ReadName').transpose()

def find(read_name, clN):
    #try:
     #   if max(m[read_name])<=0 and max(m.loc[m['ReadName'] == read_name, m.columns !='ReadName'].values)<=0:
      #      df.loc[df['ReadName'] == read_name, 'Cluster'] = 'unclustered'
       #     uncl=uncl+1
    #except:

    #print(df[(df['ReadName'] == read_name)]['Cluster'].values)



        x = m.loc[m[read_name] != -1].index
        #print(read_name)
        #print(x)
        if cl[(cl['ReadName'] == read_name)]['Cluster'].values==['NA']:
            cl.loc[cl['ReadName'] == read_name,'Cluster']=clN
        #print("childs " + x)

            for i in x:
                find(i,clN)




clN=0
uncl=0



for next_read in cl['ReadName']:
    if cl[(cl['ReadName'] == next_read)]['Cluster'].values==['NA']:
        x = m.loc[m[next_read] != -1].index
        if len(x)>0:
            clN = clN + 1
            print("new cluster"+str(clN))
            find(next_read,clN)
        else:
            uncl = uncl + 1
            print(str(next_read)+" unclustered")
            cl.loc[cl['ReadName'] == next_read, 'Cluster'] = 'unclustered'


print(str(clN)+" clusters found")
print(str(uncl)+" Reads unclassified")

#cl.to_csv("clusters_%s.csv" % edge)
#print(df)



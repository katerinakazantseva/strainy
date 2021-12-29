import csv
import pysam
from Bio import SeqIO
import pandas as pd
import numpy as np

bam="/Users/ekaterina.kazantseva/MT/test_data/test.bam"
snp="/Users/ekaterina.kazantseva/MT/test.vcf"
edge='edge_1'
R=1
#READS SNP
print ("### Reading SNPs...")
SNP_pos=[]
vcf=open(snp, "rt")
i=0
for line in vcf:

    if line.split()[0]==edge:
       i = i + 1
       SNP_pos.append(line.split()[1])

print (str(i)+" SNPs found")


#READ READS AND POSITIONS
print ("### Reading Reads...")

#SNP_pos=['9098']

#print(SNP_pos)

RName=[]
Seq=[]
Pos=[]
Start=[]
Stop=[]

bamfile = pysam.AlignmentFile(bam, "rb")
iter = bamfile.fetch(edge)



for pos in SNP_pos:
    for pileupcolumn in bamfile.pileup(edge, int(pos)-1, int(pos), stepper='samtools', min_base_quality=0, ignore_overlaps=False,
                                   ignore_orphans=False, truncate=True):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                RName.append(pileupread.alignment.query_name)
                Seq.append(pileupread.alignment.query_sequence[pileupread.query_position])
                Pos.append(pos)
                Start.append(pileupread.alignment.get_reference_positions()[0])
                Stop.append(pileupread.alignment.get_reference_positions()[len(pileupread.alignment.get_reference_positions())-1])
                #print(pileupread.alignment.query_name)
                #print(pileupread.alignment.get_reference_positions()[len(pileupread.alignment.get_reference_positions())-1])




bamfile.close()

d = {'ReadName': RName, 'Seq': Seq, 'Pos': Pos, 'Start' : Start, 'Stop' : Stop}
df = pd.DataFrame(data=d)

df = df.drop(df[(df['Seq'] == 'None')].index).reset_index()
cl=pd.DataFrame(data={'ReadName': RName})
cl=cl.drop_duplicates(subset=['ReadName'])
cl['Cluster'] = 'NA'

print (str(len(RName))+" reads found")
print(cl)



#CALCULATE DISTANCE and BUILD graph
print ("### Calculatind distances/Building adj matrix...")



print (df)
read=df.loc[ 0,'ReadName']


def distance(read1,read2,df):
    d=-1
    i=0
    for snp in SNP_pos:
        #print (snp)
        i=i+1
        try:
            b1=df.loc[( df.ReadName == read1) & (df.Pos == '%s' %snp)]['Seq'].values
            b2 =df.loc[(df.ReadName == read2) & (df.Pos == '%s' % snp)]['Seq'].values
            #print(b1)
            #print(b2)
            if len(b1)==0 and len(b2)==0:
                continue

            #elif len(b1) == 0 or len(b2) == 0:
             #   d=100
            elif len(b1)==0 or  len(b2)==0 or b1 == b2:
                if d==-1:
                    d=0
                d=d
            else:
                if d==-1:
                    d=0
                d=d+1
        except:
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



def build_adj_matrix (cl,df):
    m = pd.DataFrame(-1, index=cl['ReadName'], columns=cl['ReadName'])
    #m = pd.DataFrame(0, index=[1,2,3,4],columns=[1,2,3,4])
    for i in range(1,m.shape[1]):
        print(str(i)+"/"+str(m.shape[1])+" Reads processed \r" , end="")
        first_read=m.index[i]
        read1start=df.loc[( df.ReadName == first_read)]['Start'].values[0]
        read1stop = df.loc[(df.ReadName == first_read)]['Stop'].values[0]
        for j in range(0,i):
            second_read = m.index[j]
            read2start = df.loc[(df.ReadName == second_read)]['Start'].values[0]
            read2stop = df.loc[(df.ReadName == second_read)]['Stop'].values[0]
            #print('first '+first_read+ " "+ str(read1start)+" "+str(read1stop))
            #print('second '+second_read + " " + str(read2start) + " " + str(read2stop))
            #if read1start>read2stop or read2start>read1stop:
                #m[second_read][first_read]=-1
                #print("no intersect ")
            #else:
            m[second_read][first_read] = distance(first_read,second_read, df)

    return (m)


m=build_adj_matrix (cl,df)
#m=m.set_index('ReadName', inplace = True)
m.to_csv("adj_M_%s.csv" % edge)

#m=pd.read_csv("adj_M_edge_8.csv", index_col='ReadName')
#print(m)


print("### Removing overweighed egdes...")

def remove_edges (m, R):
    m[m >= R] = -1
    return (m)

remove_edges (m, R)

#print(m)


print("### Removing strange egdes...")

def remove_zeroedges (m):
    for read_name in cl['ReadName']:
        print(read_name)
        if max(m[read_name])<=0 and m.loc[m.index == read_name].values.max()<=0:
            m=m.drop(index=read_name)
            #modDfObj = dfObj.drop('b')
            m = m.drop(columns=read_name)
            print(str(read_name)+" removed")
            #cl.loc[cl['ReadName'] == read_name, 'Cluster'] = 'unclustered'
    return (m)

#remove_zeroedges(m)



#print()



print("### Searching connected components...")


#m=m.set_index('ReadName').transpose()

def find(read_name, clN):
    #try:
     #   if max(m[read_name])<=0 and max(m.loc[m['ReadName'] == read_name, m.columns !='ReadName'].values)<=0:
      #      df.loc[df['ReadName'] == read_name, 'Cluster'] = 'unclustered'
            #uncl=uncl+1

    #print(df[(df['ReadName'] == read_name)]['Cluster'].values)



        x = m.loc[m[read_name] != -1].index
        #print(read_name)
        #print(x)
        if cl[(cl['ReadName'] == read_name)]['Cluster'].values==['NA']:
            cl.loc[cl['ReadName'] == read_name,'Cluster']=clN
        #print("childs " + x)

            for i in x:
                find(i,clN)
    #except: pass





clN=0
uncl=0



for next_read in cl['ReadName']:
    #print(next_read)
    #try:
     #   print(1)
    #if max(m[next_read])<=0 and max(m.loc[m['ReadName'] == next_read, m.columns !='ReadName'].values)<=0:
       # print(2)
       # df.loc[df['ReadName'] == next_read, 'Cluster'] = 'unclustered'
    #except:
      #  continue
    #if next_read=='SRR13128014.1592186.1':
       #cl.loc[cl['ReadName'] == next_read, 'Cluster'] = 'unclustered'
    if cl[(cl['ReadName'] == next_read)]['Cluster'].values==['NA']:
        #print(3)
        x = m.loc[m[next_read] != -1].index
        if len(x)>0:
            clN = clN + 1
            #print("new cluster"+str(clN))
            find(next_read,clN)
        else:
            uncl = uncl + 1
            print(str(next_read)+" unclustered")
            cl.loc[cl['ReadName'] == next_read, 'Cluster'] = 'unclustered'


print(str(clN)+" clusters found in "+edge)
print(str(uncl)+" Reads unclassified")

#cl.to_csv("clusters_%s.csv" % edge)
#print(df)



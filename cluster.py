import csv
import pysam
from Bio import SeqIO
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

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

print (str(len(SNP_pos))+" SNPs found")
SNP_pos=['492','519','533','1287','1373','2746','3346','4027','4531','4597','5125','5149','5164','5239','5242','5338','5369','7232', '7383','8108', '8217']


#READ READS AND POSITIONS
print ("### Reading Reads...")





RName, Seq, Pos, Start, Stop= [],[],[],[],[]
bamfile = pysam.AlignmentFile(bam, "rb")
iter = bamfile.fetch(edge)



for pos in SNP_pos:
    for pileupcolumn in bamfile.pileup(edge, int(pos)-1, int(pos), stepper='samtools', min_base_quality=0, ignore_overlaps=False,
                                   ignore_orphans=False, truncate=True):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                RName.append(pileupread.alignment.query_name)
                Seq.append(pileupread.alignment.query_sequence[pileupread.query_position])
                #print(pileupread.alignment.query_sequence[pileupread.query_position])
                Pos.append(pos)
                Start.append(pileupread.alignment.get_reference_positions()[0])
                Stop.append(pileupread.alignment.get_reference_positions()[len(pileupread.alignment.get_reference_positions())-1])





bamfile.close()

d = {'ReadName': RName, 'Seq': Seq, 'Pos': Pos, 'Start' : Start, 'Stop' : Stop}
df = pd.DataFrame(data=d)

df = df.drop(df[(df['Seq'] == 'None')].index).reset_index()
cl=pd.DataFrame(data={'ReadName': RName})
cl=cl.drop_duplicates(subset=['ReadName'])
cl['Cluster'] = 'NA'

print (str(len(cl['ReadName']))+" reads found")
print(df)



#CALCULATE DISTANCE and BUILD graph
print ("### Calculatind distances/Building adj matrix...")



#print (df)



def distance(read1,read2,df):
    d=-1
    for snp in SNP_pos:
        try:
            b1=df.loc[( df.ReadName == read1) & (df.Pos == '%s' %snp)]['Seq'].values
            b2 =df.loc[(df.ReadName == read2) & (df.Pos == '%s' % snp)]['Seq'].values
            if b1 != b2 and b1.size and  b2.size:
                if d==-1:
                    d=0
                d=d+1
            else b1 == b2:
                if d==-1:
                    d=0
                d=d
            #else:
               # d=d+0.0001
        except:
            continue
        if d>=R:
            d=R
            break

        else:
            continue
    return (d)


def build_adj_matrix (cl,df):
    m = pd.DataFrame(-1, index=cl['ReadName'], columns=cl['ReadName'])
    for i in range(1,m.shape[1]):
        print(str(i)+"/"+str(m.shape[1])+" Reads processed \r" , end="")
        first_read=m.index[i]
        read1start=df.loc[( df.ReadName == first_read)]['Start'].values[0]
        read1stop = df.loc[(df.ReadName == first_read)]['Stop'].values[0]
        for j in range(0,i):
            second_read = m.index[j]
            read2start = df.loc[(df.ReadName == second_read)]['Start'].values[0]
            read2stop = df.loc[(df.ReadName == second_read)]['Stop'].values[0]
            #print("read1 "+ str(read1start)+"-"+str(read1stop))
            #print("read2 " + str(read2start) + "-" + str(read2stop))

            if read1start>read2stop or read2start>read1stop:
                m[second_read][first_read]=0
            #elif: добавить проверку что есть снип на пересечении
            else:
                m[second_read][first_read] = distance(first_read,second_read, df)
    return (m)


m=build_adj_matrix (cl,df)

m.to_csv("adj_M_%s.csv" % edge)

#m=pd.read_csv("adj_M_edge_8.csv",index_col='ReadName')
print(m)


print("### Removing overweighed egdes...")

def remove_edges (m, R):
    m[m >= R] = -1
    return (m)

#remove_edges (m, R)
print("### Removing srange nodes...")

nodes=[]

def change_w (m):
    m[m == 0] = 0.01
    m[m == -1] = 0
    #m[m >=R] = 0
    return (m)




print("### Searching connected components...")


def find(read_name, clN):

        x = m.loc[m[read_name].isin((0,R-1))].index
        if cl[(cl['ReadName'] == read_name)]['Cluster'].values==['NA']:
            #print("assing cl "+str(clN))
            #print(read_name)
            #print(x)
            cl.loc[cl['ReadName'] == read_name,'Cluster']=clN
            for i in x:
                find(i,clN)


clN=0
uncl=0



for next_read in cl['ReadName']:
    if cl[(cl['ReadName'] == next_read)]['Cluster'].values==['NA']:
        x = m.loc[m[next_read].isin((0,R-1))].index
        if len(x)>0:
            clN = clN + 1
            find(next_read,clN)

        else:
            uncl = uncl + 1
            print(str(next_read)+" unclustered")
            cl.loc[cl['ReadName'] == next_read, 'Cluster'] = 'unclustered'




cl.to_csv("clusters-v2_%s.csv" % edge)
#cl=pd.read_csv("clusters-v2_edge_8.csv")
#print(df)

#Calculate statistics
print("Summary for: "+edge)
print("Clusters found: "+str(clN))
print("Reads unclassified: "+str(uncl))
print("Number of reads in each cluster: ")
print(cl['Cluster'].value_counts(dropna=False))


#Build and v graph




remove_edges (m, R)
def change_w (m):
    m[m == 0] = 0.01
    m[m == -1] = 0
    #m[m >=R] = 0
    return (m)


G = nx.from_pandas_adjacency(change_w(m.transpose()))
to_remove = [(a, b) for a, b, attrs in G.edges(data=True) if attrs["weight"] >= R or attrs["weight"] == 0]
G.remove_edges_from(to_remove)

sub_graphs =list(nx.clique.find_cliques(G))
nx.draw(G, with_labels = True, width=0.1,node_color='pink', node_size=10, font_size=5)
plt.show()
#subax1 = plt.subplot(121)


#for g in nx.clique.find_cliques(G):
 #   print(len(g))
  #  if int(len(g))>60:
   #     print(g)





#print(list(nx.articulation_points(G)))
#print(list(nx.local_bridges(G,with_span=False)))
#print(list(nx.bridges(G)))
#print(nx.clustering(G))


clN=1

#for i in range(1,clN+1):
 #   print("cluster "+str(i))
  #  l=cl.loc[cl['Cluster'] == i, 'ReadName'].values
   # max_d=m.loc[m.index.isin(l), m.columns.isin(l)].max().max()
    #G = nx.from_pandas_adjacency(change_w(m.transpose().loc[m.index.isin(l), m.columns.isin(l)]))
    #print(nx.is_biconnected(G))
    #if max_d>=R:
     #   print("cluster "+str(i)+" will be splitted")
      #  to_remove = [(a,b) for a, b, attrs in G.edges(data=True) if attrs["weight"] >= R]
       # G.remove_edges_from(to_remove)
        #print(list(nx.articulation_points(G)))
        #print(list(nx.articulation_points(G)))
        #print(list(nx.local_bridges(G)))







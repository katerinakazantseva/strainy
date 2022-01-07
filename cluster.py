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


def read_snp(snp):
    SNP_pos = []
    vcf = open(snp, "rt")
    for line in vcf:
        if line.split()[0] == edge:
            SNP_pos.append(line.split()[1])

    SNP_pos = ['492', '519', '533', '1287', '1373', '2746', '3346', '4027', '4531', '4597', '5125', '5149', '5164','5239', '5242', '5338', '5369', '7232', '7383', '8108', '8217']
    print(str(len(SNP_pos)) + " SNPs found")
    return(SNP_pos)


def read_bam(file):
    bamfile = pysam.AlignmentFile(file, "rb")
    iter = bamfile.fetch(edge)
    data = {}
    for pos in SNP_pos:
        for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools', min_base_quality=0,
                                           ignore_overlaps=False,
                                           ignore_orphans=False, truncate=True):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    try:
                        data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[
                            pileupread.query_position]
                    except (KeyError):
                        data[pileupread.alignment.query_name] = {}
                        data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[
                            pileupread.query_position]

    bamfile.close()
    return(data)

def distance(read1,read2,data):
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


def build_adj_matrix (cl,data):
    m = pd.DataFrame(-1, index=cl['ReadName'], columns=cl['ReadName'])
    for i in range(1,m.shape[1]):
        print(str(i)+"/"+str(m.shape[1])+" Reads processed \r" , end="")
        first_read=m.index[i]
        for j in range(0,i):
            second_read = m.index[j]
            if len(set(data[first_read].keys()) & set(data[second_read].keys()))==0:
                m[second_read][first_read]=-1
            else:
                m[second_read][first_read] = distance(first_read,second_read, data)
    return (m)


def remove_edges (m, R):
    m[m >= R] = -1
    return (m)

def change_w (m):
    m[m == 0] = 0.01
    m[m == -1] = 0
    m[m >=R] = 0
    return (m)

#READ SNPs
print ("### Reading SNPs...")
SNP_pos=read_snp(snp)

#READ READS AND POSITIONS
print ("### Reading Reads...")
data=read_bam(bam)
cl=pd.DataFrame(data={'ReadName': data.keys()})
cl['Cluster'] = 'NA'
print (str(len(cl['ReadName']))+" reads found")


#CALCULATE DISTANCE and ADJ MATRIX
print ("### Calculatind distances/Building adj matrix...")
m=build_adj_matrix (cl,data)
m.to_csv("adj_M_%s.csv" % edge)
#m=pd.read_csv("adj_M_edge_8.csv",index_col='ReadName')
print("### Removing overweighed egdes...")
m=remove_edges (m, R)


#BUILD graph

G = nx.from_pandas_adjacency(change_w(m.transpose()))
to_remove = [(a, b) for a, b, attrs in G.edges(data=True) if  attrs["weight"] == 0]
G.remove_edges_from(to_remove)
#for g in nx.find_cliques(G):
    #print(len(g))
nx.draw(G, with_labels = False, width=0.05,node_color='pink', node_size=5, font_size=5)
#plt.show()
#plt.savefig
#plt.close()



#Calculate statistics
clN=0
uncl=0
print("Summary for: "+edge)
print("Clusters found: "+str(clN))
print("Reads unclassified: "+str(uncl))
print("Number of reads in each cluster: ")
print(cl['Cluster'].value_counts(dropna=False))












#UNUSED PART


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





#for next_read in cl['ReadName']:
 #   if cl[(cl['ReadName'] == next_read)]['Cluster'].values==['NA']:
  #      x = m.loc[m[next_read].isin((0,R-1))].index
   #     if len(x)>0:
    #        clN = clN + 1
     #       find(next_read,clN)

      #  else:
       #     uncl = uncl + 1
        #    print(str(next_read)+" unclustered")
         #   cl.loc[cl['ReadName'] == next_read, 'Cluster'] = 'unclustered'




#cl.to_csv("clusters-v2_%s.csv" % edge)
#cl=pd.read_csv("clusters-v2_edge_8.csv")
#print(df)


















from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
import gzip
import os
import networkx as nx

#from datetime import datetime
#print(datetime.today().strftime('%Y-%m-%d'))
import datetime

today = datetime.date.today()
delta = datetime.timedelta(days=60) # ~ 2 months
print(today-delta)
'''
ofile="/Users/ekaterina.kazantseva/strainy2/strainy/out_strainy/intermediate/flye_outputs/contigs.fasta"
fldir="/Users/ekaterina.kazantseva/strainy2/strainy/out_strainy/intermediate/flye_outputs/asm"       #{i}/asm/assembly.fasta"


out_file=open(ofile,"w")
clusters=[111,122,203,372,486,548,566,621,765,903]



for i in clusters:
    print(i)
    file=f"{fldir}/{i}/asm/assembly.fasta"
    for record in SeqIO.parse(file, "fasta"):
        #print(record.id)
        record.id="cluster"+str(i)
        record.name = "cluster"+str(i)
        #record.description=None
        print(record)
        #print(record.name)
        #newrecord = SeqRecord(record.seq,id=f"{i}")
        SeqIO.write(record, out_file, "fasta")
out_file.close()
'''

'''

#ofile="/Users/ekaterina.kazantseva/strainy2/strainy/out_strainy/intermediate/flye_outputs/contigs.fasta"

ov_nodes = [1,2,3,4,5]
ov_edges =[{1,2},{2,3},{3,4},{3,5}]

graph = nx.DiGraph()
graph.add_nodes_from(ov_nodes)
graph.add_edges_from(ov_edges)


for i in graph.edges:
    print(i)

cntrd_nodes=[]
for node in ov_nodes:
    suc=list(graph.successors(node))
    if len(suc)==1:
        if len(list(graph.predecessors(suc[0])))==1:
            cntrd_nodes.append((node,suc[0]))

cntrd_nodes.reverse()

for pair in cntrd_nodes:
    graph = nx.contracted_nodes(graph, pair[0], pair[1],self_loops=False)


print(graph.nodes)
'''


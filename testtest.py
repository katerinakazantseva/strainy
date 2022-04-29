edges=[1,2,3,4,5,6,7,8,9]
import multiprocessing as mp


import multiprocessing
import time


def do_something(i):
   print(edges[i])
   time.sleep(2)
   print(edges[i])

def parallel(edges):
    pool = multiprocessing.Pool(2)
    pool.map(do_something, range(0, len(edges)))
    pool.close()


import gfapy
g = gfapy.Gfa.from_file("/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/assembly_graph.gfa")

for ed in g.segments:
    if ed.name !="edge_2":
        g.rm(ed)

    if ed.name=="edge_2":
        ed.sequence=ed.sequence[1:10]

print(g.segments)




def add_child_edge(g,child):
    f = g.segments
    seq=[]
    dat = {'1':{'0': 'A', '1': 'A', '4': 'A'},'2':{'0': 'B', '1': 'B', '4': 'A'},}
    for i in f:
        if i == "edge_2":
            seq = i.sequence
            seq = list(seq)
            for key, val in dat[child].items():
                try:
                    seq[int(key)] = val

                except (ValueError):
                    continue
    seq = ''.join(seq)
    g.add_line("S\tNEW\t*")
    f = g.segments
    for i in f:
        if i == "NEW":
            new_line = i
            new_line.name = str("TEST"+str(child))
            new_line.sequence = seq
            new_line.dp = 35  # coverage
    gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/test.gfa")


add_child_edge(g,'1')
add_child_edge(g,'2')
print(g.segments)
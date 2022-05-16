


import gfapy
import multiprocessing


gfa = "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/flye_3ecoli_sim_noalt_haplo/assembly_graph.gfa"
g=gfapy.Gfa.from_file(gfa)
#g = gfapy.Gfa()

def test(j):
    #g.append("S\tNEW_%s\t*"% j)
    g.add_line("S\tNEW_%s\t*"% j)
    f = g.segments
    for i in f:
        if i == "NEW_%s" % j:
            new_line = i
            new_line.name = "test_"+str(j)
            #print(g.segment_names)




def to_neighbours(g,edge):
    to_ng=[]
    for i in g.dovetails:
        if i.from_segment.name==edge and i.from_orient=='+':
            to_ng.append(i.to_segment.name)

        if i.to_segment.name==edge and ( i.to_orient!='+'):
            to_ng.append(i.from_segment.name)

    return (to_ng)


def from_neighbours(g,edge):
    from_ng=[]
    for i in g.dovetails:
        if i.to_segment.name==edge and i.to_orient=='+':
            from_ng.append(i.from_segment.name)
            print(i)
        if i.from_segment.name == edge and (i.from_orient!='+'):
            from_ng.append(i.to_segment.name)
            print(i)


    return (from_ng)
import re
def remove_link(edge, neighbor):
    print("remove line "+str(edge)+" to "+ str(neighbor))
    edge = str(edge) + "\t"
    for i in g.edges:
        if re.search(edge, str(i)):
            if re.search(neighbor, str(i)):
                g.rm(i)

for edge in []:
    print(edge)
    to_n=to_neighbours(g,edge)
    print(to_n)
    if len(to_n)==1:
        print(from_neighbours(g,to_n[0]))
        for i in from_neighbours(g,to_n[0]):
            if len(to_neighbours(g,i))>1:
                remove_link(i, to_n[0])


print(to_neighbours(g,'edge_5'))

print(from_neighbours(g,'edge_5'))


#print(from_neighbours(g,'edge_277'))
#print(to_neighbours(g,'edge_277'))

gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/test.gfa")




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




def to_neighbours(g,edge,orient):
    to_ng=[]
    for i in g.dovetails:
        if orient=='+':
            if i.from_segment.name==edge and i.from_orient=='+':

                nei=[i.to_segment.name,i.to_orient]
                to_ng.append(nei)
        if orient == '-':
            if i.from_segment.name==edge and i.from_orient=='-':
                nei=[i.to_segment.name,i.to_orient]
                to_ng.append(nei)
    return (to_ng)



def from_neighbours(g,edge, orient):
    from_ng=[]

    for i in g.dovetails:
        if orient == '+':
            if i.to_segment.name==edge and i.to_orient=='+':
                nei = [i.from_segment.name, i.from_orient]
                from_ng.append(nei)
        if orient == '-':
            if i.to_segment.name == edge and (i.to_orient=='-'):
                nei = [i.from_segment.name, i.from_orient]
                from_ng.append(nei)
    return (from_ng)


import re
def remove_link(fr,fr_or, to, to_or):
    print("remove line: "+str(fr)+str(fr_or)+" to "+ str(to)+str(to_or))
    for i in g.dovetails:
        if i.from_segment==fr and i.from_orient==fr_or and i.to_segment==to and i.to_orient==to_or:
            g.rm(i)





def clear_links(edge):
    changed=[]
    to_n=to_neighbours(g,edge,'+')
    if len(to_n)==1:
        for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
            if len(to_neighbours(g,i[0],i[1]))>1:
                remove_link(i[0],i[1], to_n[0][0],to_n[0][1])
                changed.append(to_n[0][0])
                changed.append(i[0])
    to_n=to_neighbours(g,edge,'-')
    if len(to_n)==1:
        for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
            if len(to_neighbours(g,i[0],i[1]))>1:
                remove_link(i[0],i[1], to_n[0][0],to_n[0][1])
                changed.append(to_n[0][0])
                changed.append(i[0])
    from_n = from_neighbours(g, edge, '+')
    if len(from_n) == 1:
        for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
            if len(from_neighbours(g, i[0], i[1])) > 1:
                remove_link(from_n[0][0], from_n[0][1],i[0], i[1])
                changed.append(from_n[0][0])
                changed.append(i[0])
    from_n = from_neighbours(g, edge, '-')
    if len(from_n) == 1:
        for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
            if len(from_neighbours(g, i[0], i[1])) > 1:
                remove_link(from_n[0][0], from_n[0][1],i[0], i[1])
                changed.append(from_n[0][0])
                changed.append(i[0])
    #print(changed)
    return (changed)



def test(edge):
    changed=clear_links(edge)
    for i in changed:
        test(i)


for edge in ['edge_280']:
   test(edge)



#print(g.dovetails)
#print(from_neighbours(g,'edge_277'))
#print(to_neighbours(g,'edge_277'))

gfapy.Gfa.to_file(g,"/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/flye_3ecoli_sim_noalt/test.gfa")

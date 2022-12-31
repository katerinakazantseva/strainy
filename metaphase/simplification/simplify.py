import numpy

import gfapy
import pandas as pd
from metaphase.params import *
import re
#gaf file path

gfa ="/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/before_rr/temp2.gfa"
#gfa ="/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/before_rr/graph_before_rr.gfa"
#transformed gfa file path
gfa_transformedtest = "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/before_rr/testtest.gfa"
gfa_transformedtest1 = "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/before_rr/testtest1.gfa"
gaf_file="/Users/ekaterina.kazantseva/Downloads/aln_beforerr.gaf"
g = gfapy.Gfa.from_file(gfa)
AF=0.1
I=1000

def clear_links2(edge):
    changed=False
    to_n=to_neighbours(g,edge,'+')
    if len(to_n)==1:
        for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
            if len(to_neighbours(g,i[0],i[1]))>1:
                changed=remove_link(i[0],i[1], to_n[0][0],to_n[0][1])

    from_n = from_neighbours(g, edge, '+')
    if len(from_n) == 1:
        for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
            if len(from_neighbours(g, i[0], i[1])) > 1:
                changed=remove_link(from_n[0][0], from_n[0][1],i[0], i[1])
                #changed = True
    return (changed)



def test(g):
    #gfapy.Gfa.to_file(g, gfa_transformed)
    repeat=False
    #changed = False
    print("NEW ERA")
    #for edge in ['edge_673', 'edge_677', 'edge_806_160', 'edge_806_461', 'edge_807_113', 'edge_807_140', 'edge_807_179', 'edge_807_249']:
    for edge in g.segment_names:
        changed = clear_links2(edge)
        #print (changed)
        if changed==True:
            repeat=True

    if repeat ==True:
        test(g)



def remove_link(fr,fr_or, to, to_or):
    res=False
    for i in g.dovetails:
        if (i.from_segment == fr and i.to_segment == to) or (i.from_segment == to and i.to_segment == fr) :
            g.rm(i)
            res=True
            #print("remove line: " + str(i.from_segment) + str(i.from_orient) + " to " + str(i.to_segment) + str(i.to_orient))
    return (res)

def to_neighbours(g,edge,orient):
    to_ng=[]
    for i in g.dovetails:
            if i.from_segment.name==edge and i.from_orient=='+':
                nei=[i.to_segment.name,i.to_orient]
                to_ng.append(nei)
                #print(i)
            if i.to_segment.name==edge and i.to_orient=='-':
                nei=[i.from_segment.name,i.from_orient]
                to_ng.append(nei)
                #print(i)
    return (to_ng)



def from_neighbours(g,edge, orient):
    from_ng=[]
    for i in g.dovetails:

            if i.to_segment.name==edge and i.to_orient=='+':
                nei = [i.from_segment.name, i.from_orient]
                from_ng.append(nei)

            if i.from_segment.name == edge and (i.from_orient=='-'):
                nei = [i.to_segment.name, i.to_orient]
                from_ng.append(nei)
    return (from_ng)


def add_link(fr, fr_or, to, to_or,w):
    link = 'L	%s	%s	%s	%s	0M	ex:i:%s' % (fr, fr_or, to, to_or, w)
    #print(link)
    try:
        g.add_line(link)
        print("link added from %s %s to %s %s" % (fr, fr_or, to, to_or))
    except(gfapy.NotUniqueError): pass
        #print("dd")






unitigs=["edge_177", 'edge_159','edge_354', 'edge_395']
unitigs=g.segment_names
#unitigs=[]


#unitigs=["edge_114"]

def simple(g,unitig,remove_clusters):
    n={}
    opt=None
    #print(to_neighbours(g, unitig, '+'))
    #print(from_neighbours(g, unitig, '+'))
    N=min(len(to_neighbours(g, unitig, '+')),len(from_neighbours(g, unitig, '+')))
    if N==len(to_neighbours(g, unitig, '+')):
        #print("opt1")
        opt=1
        set1=to_neighbours(g, unitig, '+')
        set2=from_neighbours(g, unitig, '+')
    else:
        #print("opt2")
        opt = 2
        set2=to_neighbours(g, unitig, '+')
        set1=from_neighbours(g, unitig, '+')

    for i in set1:
        nei={}
        if len(i[0].split('_'))>2:
            clN=i[0].split('_')[-1]
            edge="edge_%s" % i[0].split('_')[-2]
        else:
            edge=i[0]
            clN=None
        try:
            cl=pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
        except(FileNotFoundError): break
        #print("cont")
        if clN!=None:

            if type(cl['Cluster'][0]) != str:

                reads = list(cl.loc[cl['Cluster'] == int(clN), 'ReadName'])
            else:
                reads = list(cl.loc[cl['Cluster'] == clN, 'ReadName'])
        else:
            reads = list(cl['ReadName'])
        #print(type(clN))
        #print(reads)
        #print(from_neighbours(g, unitig,'+'))
        for j in set2:
            #print(j)
            if len(j[0].split('_')) > 2:
                clN2 = j[0].split('_')[-1]
                edge2 = "edge_%s" % j[0].split('_')[-2]
            else:
                edge2 = j[0]
                clN2 = None
            try:
                cl2 = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge2, I, AF), keep_default_na=False)
            except(FileNotFoundError): break
            #print(clN2)
            if clN2 != None:
                if type(cl2['Cluster'][0])!=str:
                    reads2 = list(cl2.loc[cl2['Cluster'] == int(clN2), 'ReadName'])
                else:
                    reads2 = list(cl2.loc[cl2['Cluster'] == clN2, 'ReadName'])
            else:
                reads2 = list(cl2['ReadName'])

            if len(set(reads).intersection(set(reads2)))>0:
                nei[j[0]]=len(set(reads).intersection(set(reads2)))

            if len(reads)==0 or len(reads2)==0:
                print("check")
            #print(reads2)

        n[i[0]]=nei
    #print(n)
    for k in n.keys():
        if len(k.split('_')) == 2:
            for k2 in n.keys():
                if len(k2.split('_')) > 2 and k.split('_')[0] == k2.split('_')[0]:
                    for i in n[k2].keys():
                        try:
                            n[k][i]=n[k][i]-n[k2][i]
                        except(KeyError):pass
        else:
         for i in n[k].keys():
             if len(i.split('_')) == 2:
                 for i2 in n[k].keys():
                     if len(i2.split('_')) > 2 and i.split('_')[0] == i2.split('_')[0]:
                         try:
                            n[k][i] = n[k][i] - n[k][i2]
                         except(KeyError): pass


    #print(n)
    for i in n.keys():
        if len(n[i])==0:
            set1=[]
            set2=[]
    if len(n)==0:
        set1 = []
        set2 = []


    all_n = []
    N = 0
    for i in set1+set2:
        all_n.append(i[0])
    for i in set1:
        try:
            nei=n[i[0]]
        except(KeyError): pass

        '''
        for k in nei.keys():
            if len(k.split('_'))==2:
                for k2 in nei.keys():
                    if len(k2.split('_'))>2 and k.split('_')[0]==k2.split('_')[0]:
                        nei[k]=nei[k]-nei[k2]
        #print(nei)

        '''



        if len(nei)>0:
            j = max(nei, key=nei.get)
            N = N + 1
            #for j in nei.keys():

            #print(j)

            #print("PAIR "+ str(i)+" and "+ str(j))
            #print(N)

            #try:
            if 1==1:
                    try:
                        all_n.remove(i[0])
                        all_n.remove(j)
                    except(ValueError): pass
                    g.add_line("S\t%s_%s\t*" % (unitig, N))
                    print("CREATE S\t%s_%s\t*" % (unitig, N))
                    if opt==2:
                    #
                        add_link("%s_%s" % (unitig, N), '+', j, '+', 1)
                        remove_link(j,'+', unitig, '+')

                        add_link(i[0], '+', "%s_%s" % (unitig, N), '+', 1)
                        remove_link(unitig, '+', i[0], '+')
                    elif opt==1:
                        add_link(j, '+', "%s_%s" % (unitig, N), '+', 1)
                        add_link("%s_%s" % (unitig, N), '+', i[0], '+', 1)
                        remove_link(j, '+', unitig, '+')
                        remove_link(unitig, '+', i[0], '+')


            #except (gfapy.error.NotUniqueError):
                    #continue

    #print("set")
    #print(all_n)
    if len(all_n)>0:
        for i in all_n:
            mx=0
            N=0
            for k in n.keys():
                N=N+1
                try:
                    if n[k][i]>mx:
                        pair=k
                        mx=n[k][i]
                except: pass

            add_link("%s_%s" % (unitig, N), '+',i , '+', 1)
            #print("add "+str(i)+" "+"%s_%s" % (unitig, N))










    if N>0:
        for i in g.segments:
            if i == unitig:
                cov=i.dp
                seq=i.sequence


        for j in range(1,N+1):

            for i in g.segments:
                if i.name == "%s_%s"  % (unitig, j):
                    i.dp=round(cov/N)
                    i.sequence=seq
        remove_clusters.append(unitig)
    return (remove_clusters)



remove_clusters=[]



for unitig in unitigs:

    #print(unitig)
    if len(to_neighbours(g, unitig,'+'))>1 and len(from_neighbours(g, unitig,'+'))>1:

        print(unitig)
        remove_clusters=simple(g,unitig,remove_clusters)
        print()



for unitig in remove_clusters:
    g.rm(unitig)
    #print("REMOVE")
    print(unitig)
    for d in g.dovetails:
        if d.from_segment ==unitig or d.to_segment==unitig:
            g.rm(d)
            #print(d)
            #print("removed")

#print(g.try_get_segment("edge_444"))
#print(g.dovetails)

'''
    remove=True
    print("REMOVE")
    print(unitig)
    if len(to_neighbours(g,unitig,'+'))>0 and len(to_neighbours(g,unitig,'+'))>0:
        N=1
        while ("%s_%s" % (unitig, N) in g.segment_names) == True:
            N=N+1

        for i in g.segments:
                if i.name == unitig:
                    cov = i.dp / (N + 1)
                    i.dp=round(cov)
        for j in range(1,N+1):
            for i in g.segments:
                if i.name == "%s_%s"  % (unitig, j):
                    i.dp=round(cov)


    else:
        g.rm(unitig)'''





gfapy.Gfa.to_file(g,gfa_transformedtest)
test(g)
gfapy.Gfa.to_file(g,gfa_transformedtest1)


#print(gfapy.Gfa.try_get_line())



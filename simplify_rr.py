import numpy

import gfapy
import pandas as pd
from params import *
import re
#gaf file path

gfa ="/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/minigraph/test_split_1.gfa"
#gfa ="/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/before_rr/graph_before_rr.gfa"
#transformed gfa file path
gfa_transformedtest = "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/minigraph/test_split_rr.gfa"
gfa_transformedtest1 = "/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/minigraph/test_split_rr+sm.gfa"
gaf_file="/Users/ekaterina.kazantseva/Downloads/aln_beforerr.gaf"

AF=0.1
I=1000
''''''


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








N_new_unitigs={}

def simple(g,unitig,remove_clusters):
    n={}
    opt=None
    created=False
    N=min(len(to_neighbours(g, unitig, '+')),len(from_neighbours(g, unitig, '+')))
    if N==len(to_neighbours(g, unitig, '+')):
        opt=1
        set1=to_neighbours(g, unitig, '+')
        set2=from_neighbours(g, unitig, '+')
    else:
        opt = 2
        set2=to_neighbours(g, unitig, '+')
        set1=from_neighbours(g, unitig, '+')

    for i in set1:
        nei={}
        if len(i[0].split('_'))>1: #change for two for edges (metaflye) 1 for mini
            clN=i[0].split('_')[-1]
            #edge="edge_%s" % i[0].split('_')[-2]
            edge = i[0].split('_')[0]
        else:
            edge=i[0]
            clN=None


        try:
            cl=pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
        except(FileNotFoundError): break
        if clN!=None:

            if type(cl['Cluster'][0]) != str:

                reads = list(cl.loc[cl['Cluster'] == int(clN), 'ReadName'])
            else:
                reads = list(cl.loc[cl['Cluster'] == clN, 'ReadName'])
        else:
            reads = list(cl['ReadName'])
        #print(type(clN))
        print(reads)
        #print(from_neighbours(g, unitig,'+'))
        for j in set2:
            #print(j)
            if len(j[0].split('_')) > 1: #change for two
                clN2 = j[0].split('_')[-1]
                #edge2 = "edge_%s" % j[0].split('_')[-2]
                edge2 = j[0].split('_')[0]
            else:
                edge2 = j[0]
                clN2 = None
            try:
                cl2 = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge2, I, AF), keep_default_na=False)
            except(FileNotFoundError): break



            if clN2 != None:
                if type(cl2['Cluster'][0])!=str:
                    reads2 = list(cl2.loc[cl2['Cluster'] == int(clN2), 'ReadName'])
                    print("1")
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

    print(n)
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


    print(n)
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
        print(i)
        try:
            nei=n[i[0]]

        except(KeyError):
            nei={}

        print(nei)
        if len(nei)>0:
            print("yes")
            j = max(nei, key=nei.get)
            N = N + 1
            if 1==1:
                    try:
                        all_n.remove(i[0])
                        all_n.remove(j)
                    except(ValueError): pass
                    g.add_line("S\t%s_%s\t*" % (unitig, N))
                    print("CREATE S\t%s_%s\t*" % (unitig, N))
                    created==True
                    if opt==2:
                        add_link("%s_%s" % (unitig, N), '+', j, '+', 1)
                        remove_link(j,'+', unitig, '+')
                        add_link(i[0], '+', "%s_%s" % (unitig, N), '+', 1)
                        remove_link(unitig, '+', i[0], '+')
                    elif opt==1:
                        try:
                            n_j=N_new_unitigs[j]
                            for nj in range (1,n_j):
                                add_link("%s_%s" % (j,nj), '+', "%s_%s" % (unitig, N), '+', 1)
                        except (KeyError):
                                #print("cccc")
                                add_link(j, '+', "%s_%s" % (unitig, N), '+', 1)

                        try:
                            n_i=N_new_unitigs[i[0]]
                            for ni in range (1,n_i):
                                add_link("%s_%s" % (unitig, N), '+', "%s_%s" % (i[0], ni), '+', 1)
                        except (KeyError):
                                #print("cccc22")
                                add_link("%s_%s" % (unitig, N), '+', i[0], '+', 1)

                        remove_link(j, '+', unitig, '+')
                        remove_link(unitig, '+', i[0], '+')


            #except (gfapy.error.NotUniqueError):
                    #continue

    #print(opt)
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

    '''
    i = g.try_get_segment(unitig)
    seq = i.sequence
    ln=len(seq)
    cov = int(i.dp)

    if len(all_n)==0 and created==False and ln<15000 and len(to_neighbours(g, unitig, '+'))==len(from_neighbours(g, unitig, '+')):
        print("New code")
        for i in n:
            N=N+1
            g.add_line("S\t%s_%s\t*" % (unitig, N))
            print("CREATE S\t%s_%s\t*" % (unitig, N))
            if opt == 1:
                for d in g.dovetails:
                    if d.from_segment==unitig and d.to_segment==i:
                        #add_link("%s_%s" % (unitig, N), '+', d.to_segment.name, '+', 1)
                        link=str(d).replace(unitig,"%s_%s" % (unitig, N) )
                        g.add_line(link)
                    if d.to_segment == unitig:
                        link=str(d).replace(unitig,"%s_%s" % (unitig, N) )
                        g.add_line(link)
                        #add_link(d.from_segment.name, '+',"%s_%s" % (unitig, N), '+', 1)


            elif opt == 2:
                for d in g.dovetails:

                    if d.to_segment==unitig and d.from_segment==i:
                        #print(d)
                        #add_link(d.from_segment.name, '+',"%s_%s" % (unitig, N), '+',  1)
                        link=str(d).replace(unitig,"%s_%s" % (unitig, N) )
                        g.add_line(link)
                    if d.from_segment == unitig:
                        #print(d)
                        #add_link("%s_%s" % (unitig, N), '+',d.to_segment.name, '+', 1)
                        link=str(d).replace(unitig,"%s_%s" % (unitig, N) )
                        g.add_line(link) 
                                            '''
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

    N_new_unitigs[unitig]=N


    return (remove_clusters)



remove_clusters=[]


#unitigs=["s325_1"]

g = gfapy.Gfa.from_file(gfa, vlevel=0)
unitigs=g.segment_names
for unitig in unitigs:
    try:
        i = g.try_get_segment(unitig)
        cov = int(i.dp)
        if len(to_neighbours(g, unitig,'+'))>1 and len(from_neighbours(g, unitig,'+'))>1: #and ((len(unitig.split('_'))==2 and unitig.split('_')[-1]==1) or len(unitig.split('_'))==1):

            print(unitig)
            remove_clusters=simple(g,unitig,remove_clusters)
            print()
    except:
        pass

    gfapy.Gfa.to_file(g, gfa_transformedtest)



for unitig in remove_clusters:
    g.rm(unitig)
    #print("REMOVE")
    print(unitig)
    for d in g.dovetails:
        if d.from_segment ==unitig or d.to_segment==unitig:
            g.rm(d)
            #print(d)
            #print("removed")



from simplify_links import *
gfapy.Gfa.to_file(g,gfa_transformedtest)

simplify_links(g)
gfapy.Gfa.to_file(g,gfa_transformedtest1)





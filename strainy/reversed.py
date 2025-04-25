import gfapy
from collections import deque
#import pandas as pd




def format_links(graph):
    links={}
    for line in graph.edges:
        l1 = str(line).split()[1]
        l2 = str(line).split()[3]
        try:
            links[l1]
        except KeyError:
            links[l1]={}
            links[l1]["+"] = []
            links[l1]["-"] = []
        try:
            links[l2]
        except KeyError:
            links[l2]={}
            links[l2]["+"] = []
            links[l2]["-"] = []

        if str(line).split()[2] == "+" and str(line).split()[4] == "+":
            try:
                l1_n=links[l1]["+"]
                if l2 not in l1_n:
                    l1_n.append(l2)
                    links[l1]["+"]=l1_n
            except KeyError:
                print("check")
            try:
                l2_n=links[l2]["+"]
                if l1 not in l2_n:
                    l2_n.append(l1)
                    links[l2]["+"]=l2_n
            except KeyError:
                print("check")
        else:
            try:
                l1_n=links[l1]["-"]
                if l2 not in l1_n:
                    l1_n.append(l2)
                    links[l1]["-"]=l1_n
            except KeyError:
                print("check")
            try:
                l2_n=links[l2]["-"]
                if l1 not in l2_n:
                    l2_n.append(l1)
                    links[l2]["-"]=l2_n
            except KeyError:
                print("check")
    return links




def bfs(links, edges,unknown):
    pos = []
    neg = []
    #unknown = []
    visited = []
    s = edges[0]
    pos.append(s)
    # Create a queue for BFS
    q = deque()

    # Initially mark all the vertices as not visited
    # When we push a vertex into the q, we mark it as
    # visited

    visited.append(s)
    # Mark the source node as visited and enqueue it
    q.append(s)

    # Iterate over the queue
    while q:
        curr = q.popleft()
        nei_for=links[curr]["+"]
        nei_rev = links[curr]["-"]
        nei=nei_rev+nei_for
        for x in nei:
            if (x in nei_for and curr in pos) or (x in nei_rev and curr in neg) :
                if x in neg and x not in unknown:
                    unknown.append(curr)
                    pos.append(curr)
                    neg.append(curr)
                elif x not in pos : #
                    pos.append(x)
            if (x in nei_for and curr in neg) or (x in nei_rev and curr in pos) :
                if x in pos and x not in unknown:
                    unknown.append(curr)
                    pos.append(curr)
                    neg.append(curr)
                elif x not in neg : #
                    neg.append(x)

            if x not in  visited:
                visited.append(x)
                q.append(x)
    dct={}
    dct["pos"]=set(pos)
    dct["neg"] = set(neg)
    dct["unk"] = set(unknown)
    return dct





def minimize(graph):
    edges = graph.segment_names
    links = format_links(graph)
    visited = []
    unknown = []
    first_res=bfs(links, edges, [])
    first = first_res["unk"]
    candidates=[]
    for i in first:
        candidates.append(i)
        nei_for = links[i]["+"]
        nei_rev = links[i]["-"]
        nei = nei_rev + nei_for
        for n in nei:
            candidates.append(n)
    candidates=set(candidates)
    min=10000
    if len(candidates)==0:
        return first_res
    for i in candidates:
        curl=len(bfs(links, edges,[i])["neg"])
        if curl<min:
            min=curl
            res=i

    resd=bfs(links, edges,[res])
    return resd


'''
#file="/Users/ekaterina.kazantseva/strainy2/strainy/test_set/toy.gfa"
#file="/Users/ekaterina.kazantseva/Documents/MANUAL/strainy_ecoli_example/ecoli_5strain_metaflye_hap.gfa" #ok
#file="/Users/ekaterina.kazantseva/Downloads/gtest.gfa"
#file="/Users/ekaterina.kazantseva/Downloads/g1.gfa"
#file="/Users/ekaterina.kazantseva/Downloads/gtest2.gfa
#input_graph = gfapy.Gfa.from_file(file)

resd=minimize(input_graph)
print(resd["neg"])
pos = resd["pos"]
neg =  resd["neg"]
unknown=resd["unk"]
edges=input_graph.segment_names
colr=[]
for i in edges:
    if i in pos and i not in unknown:
        colr.append("red")
    elif i in neg and i not in unknown:
        colr.append("green")
    elif i in unknown:
        colr.append("black")


col = pd.DataFrame(
    {'Name': edges,
     'color': colr})

col.to_csv("/Users/ekaterina.kazantseva/Downloads/g11.csv",sep=',', index=False)
'''
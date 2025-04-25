import gfapy
import pandas as pd
file="/Users/ekaterina.kazantseva/strainy2/strainy/test_set/toy.gfa"
#file="/Users/ekaterina.kazantseva/Documents/MANUAL/strainy_ecoli_example/ecoli_5strain_metaflye_hap.gfa" #ok
file="/Users/ekaterina.kazantseva/Downloads/gtest.gfa" #ok
#file="/Users/ekaterina.kazantseva/Downloads/gtest2.gfa" #не ок
input_graph = gfapy.Gfa.from_file(file)

pos=[]
neg=[]
'''
print("1")
for line in input_graph.edges:
    l1 = str(line).split()[1]
    l2 = str(line).split()[3]


    if str(line).split()[2]=="+" and str(line).split()[4]=="+":
        if l1 not in pos and l1 not in neg:
            pos.append(l1)
        if l2 not in pos and l1 not in neg:
            pos.append(l2)

    if str(line).split()[2] == "+" and str(line).split()[4] == "-" or (str(line).split()[2] == "-" and str(line).split()[4] == "+"):
        if l1 in pos and l2 not in pos and l2 not in neg:
            neg.append(l2)
        if l2 in pos and l1 not in pos and l1 not in neg:
            neg.append(l1)

print(pos)
print(neg)

print("2")
pos=[]
neg=[]
visited=[]
edges=input_graph.segment_names
first=edges[0]
pos.append(first)
#visited.append(first)
Q=[]



'''
links={}
for line in input_graph.edges:
    l1 = str(line).split()[1]
    l2 = str(line).split()[3]
    #здесь создать записи
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


print(links)

'''
Q=edges

while Q:
    edge = Q[0]
    Q.pop(0)
    if edge in pos or edge in neg and edge not in visited:
        pos_nei=links[edge]["+"]
        for i in pos_nei:
            if i not in pos and i not in neg:
                if edge in pos:
                    pos.append(i)
                else:
                    neg.append(i)
                if i not in Q and i not in visited:
                    Q.append(i)
        neg_nei=links[edge]["-"]
        for i in neg_nei:
            if i not in pos and i not in neg:
                if edge in pos:
                    neg.append(i)
                else:
                    pos.append(i)
            #visited.append(i)
                if i not in Q and i not in visited:
                    Q.append(i)
        visited.append(edge)
    else:
        Q.append(edge)


print(pos)
print(neg)

#coord = {k: [v["Start"], v["End"]] for k, v in value.items()}
#row = pd.DataFrame({'ReadName': [key], 'Cluster': ['NA'], 'Coordinates': [coord]})
#cl = pd.concat([cl, row])

edges=input_graph.segment_names
colr=[]
for i in edges:
    if i in pos:
        colr.append("red")
    elif i in pos:
        colr.append("green")
    else:
        colr.append("black")
print(colr)
#col = pd.DataFrame(columns=['Node', 'Colour'])
col = pd.DataFrame(
    {'Name': edges,
     'color': colr})
#row = pd.DataFrame({'Edge': [edges], 'Colour': [colr]})
#col = pd.concat([col, row])
#print(row)
#print (col)
#col.to_csv("/Users/ekaterina.kazantseva/Downloads/col.csv",sep=',', index=False,)

'''

pos=[]
neg=[]
unknown=[]
visited=[]
edges=input_graph.segment_names
Q=edges
unknown=[]


max=0
for k,v in links.items():
    if len(v["+"])>max:
        max=len(v["+"])
        maxedge=k



pos.append(maxedge)

def add_pos(edge):
    for i in links[edge]["+"]:
        if i not in pos:
            pos.append(i)
            add_pos(i)

add_pos(maxedge)

for i in edges:
    if i not in pos:
        unknown.append(i)

'''
filled=0
while Q:
    edge = Q[0]
    Q.pop(0)
    if len(links[edge]["-"])==0 and filled==0:
        filled=1
        if edge not in pos or edge  not in neg:
            pos.append(edge)
            for i in links[edge]["+"]:
                pos.append(i)
    else:
        unknown.append(edge)

'''




def check_links4(edge,pos,neg,unknown):
  a="NA"
  b="NA"
  if len([ i for i in  links[edge]["+"] if i in pos])==len(links[edge]["+"]) and len(links[edge]["+"])!=0 :
      a="POS"
  elif len([ i for i in  links[edge]["+"] if i in neg])==len(links[edge]["+"]) and len(links[edge]["+"])!=0:
      a="NEG"
  elif len([i for i in links[edge]["+"] if i in unknown]) == len(links[edge]["+"]) and len(links[edge]["+"])!=0:
    a = "unk"
  if len([i for i in links[edge]["-"] if i in pos]) == len(links[edge]["-"])  and len(links[edge]["-"])!=0:
    b = "POS"
  elif len([i for i in links[edge]["-"] if i in neg]) == len(links[edge]["-"]) and len(links[edge]["-"])!=0:
    b = "NEG"
  elif len([i for i in links[edge]["-"] if i in unknown]) == len(links[edge]["-"]) and len(links[edge]["-"])!=0:
    b = "unk"
  return (a,b)

def check_links(edge,pos,neg,unknown):
  if len(links[edge]["+"])==0:
      a = "NA"
  elif len([ i for i in  links[edge]["+"] if i in pos])==len([ i for i in  links[edge]["+"] if i in neg]):
      a= "unk"
  elif len([i for i in links[edge]["+"] if i in pos]) > len([i for i in links[edge]["+"] if i in neg]):
      a = "POS"
  else:
      a="NEG"
  if len(links[edge]["-"])==0:
      b = "NA"
  elif len([ i for i in  links[edge]["-"] if i in pos])==len([ i for i in  links[edge]["-"] if i in neg]):
      b= "unk"
  elif len([i for i in links[edge]["-"] if i in pos]) > len([i for i in links[edge]["-"] if i in neg]):
      b = "POS"
  else:
      b="NEG"

  return (a,b)


visited=[]
contr=[]

def update(edge):
    if edge not in visited and edge in unknown:
        a, b=check_links(edge,pos,neg, unknown)
        print(edge)
        print(links[edge])
        print(a,b)
        nei=[]
        if  (a=="POS" and b=="NA") : #(a=="POS" and b=="unk") or
            pos.append(edge)
            unknown.remove(edge)
            nei = links[edge]["-"]
            visited.append(edge)
            print("move to pos")
            print()
            for i in nei:
                update(i)

        elif  (a=="NEG" and b=="NA") : #(a=="NEG" and b=="unk") or
            neg.append(edge)
            try:
                unknown.remove(edge)
            except ValueError:
                pass
            nei = links[edge]["+"]
            visited.append(edge)
            print("move to neg")
            print()
            for i in nei:
                update(i)

        elif (b=="POS" and a=="NA") : #(a=="unk" and b=="POS") or
            neg.append(edge)
            unknown.remove(edge)
            nei = links[edge]["-"]
            visited.append(edge)
            print("move to neg")
            print()
            for i in nei:
                update(i)

        elif (a=="POS" and b=="NEG"): #or (a=="POS" and b=="POS"):
            pos.append(edge)
            unknown.remove(edge)
            nei = links[edge]["-"]
            visited.append(edge)
            print("move to pos")
            print()
            for i in nei:
                update(i)
        elif (a=="unk" and b=="NEG"): #or (a=="POS" and b=="POS"):
            pos.append(edge)
            unknown.remove(edge)
            nei = links[edge]["-"]
            visited.append(edge)
            print("move to pos")
            print()
            for i in nei:
                update(i)
        print()


print(pos)
print(neg)
print(unknown)

q=unknown.copy()
while q:
    print("next")
    edge=q[0]
    print(edge)
    update(edge)
    q.pop(0)




print(pos)
print(neg)
print(unknown)

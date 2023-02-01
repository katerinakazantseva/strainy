import gfapy
from strainy.params import *
import logging


logger = logging.getLogger()


def to_neighbours(g,edge,orient):
    if orient =="+":
        to_ng=[]
        for i in g.dovetails:
                if i.from_segment.name==edge and i.from_orient=='+':
                    nei=[i.to_segment.name,i.to_orient]
                    to_ng.append(nei)
                if i.to_segment.name==edge and i.to_orient=='-':
                    if i.from_orient == "+":
                        nei=[i.from_segment.name,"-"]
                    elif i.from_orient == "-":
                        nei=[i.from_segment.name,"+"]
                    to_ng.append(nei)

    else:
        to_ng = []
        for i in g.dovetails:
            if i.to_segment.name == edge and i.to_orient == '+':
                if i.from_orient=="+":
                    nei = [i.from_segment.name, "-"]
                elif i.from_orient=="-":
                    nei = [i.from_segment.name, "+"]
                to_ng.append(nei)

            if i.from_segment.name == edge and (i.from_orient == '-'):
                nei = [i.to_segment.name, i.to_orient]
                to_ng.append(nei)
    return (to_ng)



def from_neighbours(g,edge, orient):
    if orient == '+':
        from_ng=[]
        for i in g.dovetails:

                if i.to_segment.name==edge and i.to_orient=='+':
                    nei = [i.from_segment.name, i.from_orient]
                    from_ng.append(nei)

                if i.from_segment.name == edge and (i.from_orient=='-'):
                    if i.to_orient == "+":
                        nei = [i.to_segment.name, "-"]
                    elif i.to_orient == "-":
                        nei = [i.to_segment.name, "+"]
                    from_ng.append(nei)

    else:
        from_ng = []
        for i in g.dovetails:
            if i.from_segment.name == edge and i.from_orient == '+':
                if i.to_orient == "+":
                    nei = [i.to_segment.name, "-"]
                elif i.to_orient == "-":
                    nei = [i.to_segment.name, "+"]
                from_ng.append(nei)
            if i.to_segment.name == edge and i.to_orient == '-':
                nei = [i.from_segment.name, "-"]
                from_ng.append(nei)
    return (from_ng)

def clear_links(edge,g):
    try:
        i=g.try_get_segment(edge)
        edge_cov = int(i.dp)
        if edge_cov==0:
            edge_cov=0.1
    except:
        edge_cov = 0.1
        pass
    changed=False
    to_n=to_neighbours(g,edge,'+')
    try:
        if len(to_n)==1 and g.try_get_segment(to_n[0][0]).dp/edge_cov<cov_ratio:
            for i in from_neighbours(g,to_n[0][0],to_n[0][1]):
                if i[0] != edge:
                    if len(to_neighbours(g,i[0],i[1]))>1:
                        changed=remove_link(i[0],i[1], to_n[0][0],to_n[0][1],g)
    except:
        pass

    from_n = from_neighbours(g, edge, '+')
    try:
        if len(from_n) == 1 and g.try_get_segment(from_n[0][0]).dp/edge_cov<cov_ratio:
            for i in to_neighbours(g, from_n[0][0], from_n[0][1]):
                if i[0]!=edge:
                    if len(from_neighbours(g, i[0], i[1])) > 1:
                        changed=remove_link(from_n[0][0], from_n[0][1],i[0], i[1],g)
    except:
        pass
    return (changed)

def remove_link(fr,fr_or, to, to_or,g):
    res=False
    for i in g.dovetails:
        if (i.from_segment == fr and i.to_segment == to) or (i.from_segment == to and i.to_segment == fr):
            g.rm(i)
            res=True
            logger.debug("remove line: " + str(i.from_segment) + str(i.from_orient) + " to " + str(i.to_segment) + str(i.to_orient))
    return (res)



def remove_zero_cov(g):
    to_remove=[]
    for unitig in g.segment_names:
        i = g.try_get_segment(unitig)
        try:
            cov = i.dp
        except:
            pass
        if cov==0:
            g.rm(unitig)
            to_remove.append(unitig)
            logger.debug("REMOVE " + unitig)
    for d in g.dovetails:
        if d.from_segment in to_remove or d.to_segment in to_remove:
            g.rm(d)
    return(g)

def simplify_links(g):
    if minigraph==True:
        g=remove_zero_cov(g)
    repeat=False
    for edge in g.segment_names:
        changed = clear_links(edge,g)
        if changed==True:
            repeat=True
    if repeat ==True:
        simplify_links(g)





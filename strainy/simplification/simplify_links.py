import gfapy
from strainy.params import *
import logging


logger = logging.getLogger()


def to_neighbours(g,edge,orient):
    to_ng = []
    segment = g.try_get_segment(edge)

    if orient == "+":
        for link in segment.dovetails_R:
            if link.from_segment.name == edge and link.from_orient == "+":
                nei = [link.to_segment.name, link.to_orient]
                to_ng.append(nei)
            if link.to_segment.name == edge and link.to_orient=="-":
                if link.from_orient == "+":
                    nei=[link.from_segment.name,"-"]
                elif link.from_orient == "-":
                    nei=[link.from_segment.name,"+"]
                to_ng.append(nei)

    else:
        for link in segment.dovetails_L:
            if link.to_segment.name == edge and link.to_orient == "+":
                if link.from_orient=="+":
                    nei = [link.from_segment.name, "-"]
                elif link.from_orient=="-":
                    nei = [link.from_segment.name, "+"]
                to_ng.append(nei)

            if link.from_segment.name == edge and (link.from_orient == "-"):
                nei = [link.to_segment.name, link.to_orient]
                to_ng.append(nei)
    return to_ng


def from_neighbours(g,edge, orient):
    from_ng = []
    segment = g.try_get_segment(edge)

    if orient == "+":
        for link in segment.dovetails_L:
            if link.to_segment.name == edge and link.to_orient == "+":
                nei = [link.from_segment.name, link.from_orient]
                from_ng.append(nei)

            if link.from_segment.name == edge and (link.from_orient == "-"):
                if link.to_orient == "+":
                    nei = [link.to_segment.name, "-"]
                elif link.to_orient == "-":
                    nei = [link.to_segment.name, "+"]
                from_ng.append(nei)

    else:
        for link in segment.dovetails_R:
            if link.from_segment.name == edge and link.from_orient == "+":
                if link.to_orient == "+":
                    nei = [link.to_segment.name, "-"]
                elif link.to_orient == "-":
                    nei = [link.to_segment.name, "+"]
                from_ng.append(nei)
            if link.to_segment.name == edge and link.to_orient == "-":
                nei = [link.from_segment.name, "-"]
                from_ng.append(nei)

    return from_ng


def clear_links(edge,g):
    edge_cov = g.try_get_segment(edge).dp
    if edge_cov == 0:
        edge_cov = 0.1

    changed = False
    to_n = to_neighbours(g, edge, "+")
    try:
        if len(to_n) == 1 and g.try_get_segment(to_n[0][0]).dp / edge_cov < cov_ratio:
            siblings = [u for u in from_neighbours(g,to_n[0][0], to_n[0][1]) if u[0] != edge]
            if len(siblings) == 1:
                for i in siblings:
                    if len(to_neighbours(g,i[0],i[1])) > 1:
                        changed = remove_link(i[0],i[1], to_n[0][0],to_n[0][1],g)
    except TypeError as e:
        logger.error(f'Caught {e} for {g.try_get_segment(to_n[0][0])}')
    try:
        from_n = from_neighbours(g, edge, "+")
        if len(from_n) == 1 and g.try_get_segment(from_n[0][0]).dp / edge_cov < cov_ratio:
            siblings = [u for u in to_neighbours(g, from_n[0][0], from_n[0][1]) if u[0] != edge]
            if len(siblings) == 1:
                for i in siblings:
                    if len(from_neighbours(g, i[0], i[1])) > 1:
                        changed = remove_link(from_n[0][0], from_n[0][1],i[0], i[1],g)
    except TypeError as e: 
        logger.error(f'Caught {e} for {g.try_get_segment(from_n[0][0])}')
    return changed


def remove_link(fr, fr_or, to, to_or, g):
    if fr == to:    #Weird self loops that should not happen, it seems it breaks somwthing in gfapy
        return False

    res = False
    first_seg = g.try_get_segment(fr)

    for link in first_seg.dovetails_L + first_seg.dovetails_R:
        if (link.from_segment == fr and link.to_segment == to) or (link.from_segment == to and link.to_segment == fr):
            g.rm(link)
            res = True
            logger.debug(f"remove line: {link.from_segment}{link.from_orient} to {link.to_segment}{link.to_orient}")

    return res


def remove_zero_cov(g):
    to_remove = []
    for unitig in g.segment_names:
        cov = g.try_get_segment(unitig).dp
        if cov == 0:
            g.rm(unitig)
            to_remove.append(unitig)
            logger.debug("REMOVE " + unitig)
    for d in g.dovetails:
        if d.from_segment in to_remove or d.to_segment in to_remove:
            g.rm(d)
    return(g)


def simplify_links(g):
    logger.debug("Simplification iteration")
    if minigraph == True:
        g = remove_zero_cov(g)
    repeat = False
    for edge in g.segment_names:
        changed = clear_links(edge,g)
        if changed == True:
            repeat = True
    if repeat == True:
        simplify_links(g)

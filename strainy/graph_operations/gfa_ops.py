import gfapy
import networkx as nx
import logging
import re
from strainy.logging import set_thread_logging

"""
This contains functions for operation with graph of gfa format:
1. add_link: Adds a link between specified segments in the graph
2. add_edge: Adds an empty(no sequence) segment with the specified name and coverage to the graph
3. gfa_to_nx: Сonverts the graph from the gfa format to nx (networkx) format
4. from_pandas_adjacency_notinplace: Workaround for networkx.from_pandas_adjacency issue https://github.com/networkx/networkx/issues/7407
5. clean_graph: Cleans graph from selflinks, and add "A" sequence to 0-length edges
"""
logger = logging.getLogger()




def add_link(graph, fr, fr_or, to, to_or, w):
    """
    Adds a link between specified segments in the graph
    Parameters:
        graph (gfa): graph
        fr, to (string): names of segments to be linked (from and to)
        fr_or,fr_or (string): orientation of segments to be linked (from and to)
        w: weight of the link
    """
    #check if segments exist before connecting
    if graph.segment(fr) is None or graph.segment(to) is None:
        return
    link = f"L\t{fr}\t{fr_or}\t{to}\t{to_or}\t0M\tex:i:{w}"
    try:
        graph.add_line(link)
        logger.debug("link added: " + link)
    except(gfapy.NotUniqueError):   #link already exists
        pass


def add_edge(graph,edge, clN, cov):
    #TODO remove edge,clN from parameters, use name instead
    #TODO add seq
    """
    Adds an empty(no sequence) segment with the specified name and coverage to the graph
    Parameters:
        graph (gfa): graph
        name (string): name of the edge to be created
        cov (coverage): coverage of the edge to be created
    Returns:
        gfa edge
    """
    graph.add_line("S\t%s_%s\t*" % (edge, clN))
    new_line = graph.try_get_segment("%s_%s" % (edge, clN))
    new_line.name = str(edge) + "_" + str(clN)
    new_line.sid = str(edge) + "_" + str(clN)
    new_line.dp = cov  # TODO: what to do with coverage?
    return new_line



def gfa_to_nx(g):
    """
    Сonverts the graph from the gfa format to nx (networkx) format
    Parameters:
        graph (gfa): gfa graph
    Returns:
        graph (nx): nx graph
    """
    G = nx.Graph()
    for i in g.segment_names:
        G.add_node(i)
    for i in g.dovetails:
        G.add_edge(i.from_segment.name, i.to_segment.name)
    return G


def from_pandas_adjacency_notinplace(df, create_using=None):
    """
    This function is exactly the same as networkx.from_pandas_adjacency,
    with the exception of 'copy=True' argument in relabel_nodes.
    This is because copy=False (default option) implies that the graph
    can be relabeled in place, which is not always possible
    https://github.com/networkx/networkx/issues/7407
    """
    try:
        df = df[df.index]
    except Exception as err:
        missing = list(set(df.index).difference(set(df.columns)))
        msg = f"{missing} not in columns"
        raise nx.NetworkXError("Columns must match Indices.", msg) from err

    A = df.values
    G = nx.from_numpy_array(A, create_using=create_using)

    G = nx.relabel.relabel_nodes(G, dict(enumerate(df.columns)), copy=True)
    return G


def clean_graph(g):
    """
    Сlears gfa graph - deletes edges with zero length, and self links
    Parameters:
        graph (gfa): gfa graph
    Returns:
        graph (gfa): gfa graph
    """
    for line in g.dovetails:
        if line.from_segment == line.to_segment: #TODO do not self links
            g.rm(line)
        #if g.segment(line.from_segment).virtual == True or g.segment(line.to_segment).virtual == True:
        #    g.rm(line)
    for seq in g.segments:  #TODO do not create o len unitigs
        if len(seq.sequence) == 0:
            seq.sequence = "A"
            #seq.dp = 0
    for path in g.paths:
        g.rm(path)
    return g




import gfapy
import networkx as nx
import logging

from strainy.logging import set_thread_logging


logger = logging.getLogger()

def add_link(graph, fr, fr_or, to, to_or, w):
    """
     Add gfa links between unitigs
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


def gfa_to_nx(g):
    G = nx.Graph()
    for i in g.segment_names:
        G.add_node(i)
    for i in g.dovetails:
        G.add_edge(i.from_segment.name, i.to_segment.name)
    return(G)


def from_pandas_adjacency_notinplace(df, create_using=None):
    """
    This function is exactly the same as networkx.from_pandas_adjacency,
    with the exception of 'copy=True' argument in relabel_nodes.
    This is because copy=False (default option) implies that the graph
    can be relabeled in place, which is not always possible
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

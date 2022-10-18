from karateclub import LabelPropagation
from networkx.algorithms import community
import networkx as nx


def find_communities(G):
    to_remove = [(a, b) for a, b, attrs in G.edges(data=True) if attrs["weight"] == 0]
    G.remove_edges_from(to_remove)
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    model = LabelPropagation()
    model.fit(G)
    cluster_membership = model.get_memberships()
    return (cluster_membership)
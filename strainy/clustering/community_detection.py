from karateclub import LabelPropagation
from networkx.algorithms import community
import networkx as nx


def find_communities(G):
    LabelPropagation()
    model = LabelPropagation()
    model.fit(G)
    cluster_membership = model.get_memberships()
    return cluster_membership
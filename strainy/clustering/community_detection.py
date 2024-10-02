from karateclub import LabelPropagation
from networkx.algorithms import community
import networkx as nx


def find_communities(G):
    """
        Identifies communities (clusters) within a graph using the Label Propagation algorithm.
        This function applies the Label Propagation algorithm to detect communities within the provided graph `G`.
        Each node in the graph is assigned a community label based on the algorithm, which iteratively propagates
        labels through the network until convergence.
        Returns:
            dict: A dictionary where keys are node identifiers, and values are the assigned community labels,
                  indicating the community membership of each node.
    """
    LabelPropagation()
    model = LabelPropagation()
    model.fit(G)
    cluster_membership = model.get_memberships()
    return cluster_membership
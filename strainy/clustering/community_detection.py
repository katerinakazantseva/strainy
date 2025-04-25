from karateclub import LabelPropagation
from networkx.algorithms import community
import networkx as nx
#from cdlib import algorithms #todo temove


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


    #cluster_membership = algorithms.surprise_communities(G) не работает
    #cluster_membership = algorithms.leiden(G).to_node_community_map()  #surprise
    #cluster_membership=algorithms.walktrap(G).to_node_community_map() #walk
    #cluster_membership=dict(cluster_membership)
    #for key in cluster_membership.keys():
        #cluster_membership[key] = cluster_membership[key][0]
    #print(cluster_membership)



    return cluster_membership
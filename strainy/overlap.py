from strainy.graph_operations import overlap_graph_ops
from strainy.params import StRainyArgs, init_global_args_storage
import networkx as nx




class StrainyOverlap2:
    def __init__(self,ov_nodes,ov_edges):
        self.nodes = ov_nodes
        self.edges = ov_edges
        graph = nx.DiGraph()
        graph.add_nodes_from(self.nodes)
        graph.add_edges_from(self.edges)
        self.graph = graph
    def remove_transitive(self): #move function here
        Gup=overlap_graph_ops.remove_transitive(self.graph)
        self.graph=Gup
    def merge_linear(self): #move function here
        Gup = overlap_graph_ops.remove_transitive(self.graph)
        self.graph = Gup
    def vis(self):
        G=nx.nx_agraph.to_agraph(self.graph)
        G.layout(prog = "dot")
        G.draw(f"{StRainyArgs().output_intermediate}/graphs/linear_phase2.png")
    def merge_unbranching(self):
        graph=self.graph
        cntrd_nodes = []
        for node in graph.nodes:
            suc = list(graph.successors(node))
            if len(suc) == 1:
                if len(list(graph.predecessors(suc[0]))) == 1:
                    cntrd_nodes.append((node, suc[0]))
        cntrd_nodes.reverse()
        for pair in cntrd_nodes:
            graph = nx.contracted_nodes(graph, pair[0], pair[1], self_loops=False)
        self.graph=graph
        self.nodes = graph.nodes
        self.edges = graph.edges
        return cntrd_nodes



# This is a library to that contains simple functions to iterate over the graphs
# to check for their consistency, see what are the links, etc.
# I tend to work with NetworkX graphs a lot, it makes sense to bring all the functions under one module

import networkx as nx

class GraphManipulator(object):
	def __init__(self, graph):
		self.G = graph

	def get_compvals_first_neighbors(self, node, attr):
		'''Returns the values of given attribute of first neighbors of given node, assuming the attribute exists and is instantiated
		e.g. get_compvals_first_neighbors("drug_1241", "title")'''
		titles_dict = {}
		if node in self.G:
			neighbors = self.G[node]
			for neighbor in neighbors:
				titles_dict[neighbor] = self.G.node[neighbor][attr] if attr in self.G.node[neighbor] else ""
		return titles_dict

	def set_compvals_first_neighbors(self, node, attr):
		pass

	def remove_edge_condition(self, condition):
		'''Removes edges from the graph on the basis of a condition on the edge ID (Which can be a lambda function)
		e.g. edge[0][0:4] == "gene"'''
		for edge in self.G.edges():
			if condition(edge):
				self.G.remove_edge(edge[0], edge[1])

	def remove_node_condition(self, condition):
		'''Removes nodes from the graph on the basis of a condition on the node ID (Which can be a lambda function)
		e.g. node[0][0:4] == "gene"'''
		for node in self.G.nodes():
			if condition(node):
				self.G.remove_node(node)

	def remove_edge_attr_condition(self, condition, attr):
		for edge in self.G.edges():
			if condition(self.G[edge[0]][edge[1]][attr]):
				self.G.remove_edge(edge[0], edge[1])

	def remove_node_attr_condition(self, condition, attr):
		for node in self.G.nodes():
			if condition(self.G.node[node][attr]):
				self.G.remove_node(node)

	def remove_isolates(self):
		self.G.remove_nodes_from(nx.isolates(self.G))

	def get_modified_graph(self):
		return self.G


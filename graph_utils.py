# This is a library to that contains simple functions to iterate over the graphs
# to check for their consistency, see what are the links, etc.
# I tend to work with NetworkX graphs a lot, it makes sense to bring all the functions under one module

import networkx as nx
from rdflib import Graph, URIRef, Literal, Namespace, BNode
from rdflib.namespace import RDF, RDFS, DC
from copy import deepcopy
import sys, json

class GraphManipulator(object):
	''' This is a manipulator, who takes the knowledge graph (knowledge of everything that might be going on)
	and further enquires extra details from this graph given certain conditions, and then precisely manipulates certain aspects of this knowledge
	before passing it out (to subsequent recepients). Oh the manipulator, is very cunning! This is a class of such manipulators'''

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

	def get_transitive_edges(self):
		transitive_edges = {}
		for edge in self.G.edges():
			if not (edge[0], edge[1]) in transitive_edges:
				if not (edge[1], edge[0]) in transitive_edges:
					transitive_edges[edge] = False
				else:
					transitive_edges[(edge[1], edge[0])] = True
		trc = deepcopy(transitive_edges)
		for k in trc:
			if not trc[k]:
				del transitive_edges[k]
		return transitive_edges

	def condition_attrs(self, node, attr, condition=None, attr_type="set"):
		pr_terms = []
		if attr_type == "set" or attr_type == "list":
			for term in self.G.node[node][attr]:
				pr_term = condition(term) if condition else term
				pr_terms.append(pr_term)
		return pr_terms

	def get_nbrs_type(self, node, nbr_type, nbr_value):
		'''Given a node, get neighbors of a particular type/value of attr'''
		nbrs = []
		for nbr in self.G[node]:
			if self.G.node[nbr][nbr_type] == nbr_value:
				nbrs.append(nbr)
		return nbrs 

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

class GraphPresenter(object):
	''' The presenter is somewhat innocent -- he has no curiosity like the manipulator to enquire certain aspects of the graphs
	He takes the knowledge graph, or parts of it that are provided by the manipulator, as it is and prints or converts it for subsequent users.
	The manipulator uses the presenter, how sad! This is a class of such presenters.'''

	def __init__(self, graph, entities=None):
		self.G = nx.read_gpickle(graph)
		self.entities = entities

	def print_node_attr(self, nodeList, attr):
		'''Prints the attr values of each node in a nodelist in the format:
		nodeId, attrValue'''
		for node in nodeList:
			print node, self.G.node[node][attr]

	def gen_rdfG(self, namespace):
		''' Takes the Networkx Graph and generates an RDF Graph, annotated using namespaces'''
		self.rdfG = Graph()
		self.nspace = Namespace(namespace)
		def push_triples(triples):
			for triple in triples:
				self.rdfG.add(triple)
		for node in self.G.nodes():
			title_f = lambda x: Literal(x)
			type_f = lambda x: self.nspace[x.title()[0:len(x)-1]]
			push_triples(self.gen_node_triples(node, self.G.node[node], "title", DC.title, "set", title_f))
			push_triples(self.gen_node_triples(node, self.G.node[node], "type", RDF.type, "single", type_f))
		for edge in self.G.edges():
			edge_triples = self.gen_edge_triples(edge, self.G[edge[0]][edge[1]])
			#print edge_triples
			push_triples(edge_triples)
		return self.rdfG

	def gen_node_triples(self, nodeId, nodeInfo, attr="title", pred=DC.title, attr_type="set", condition=None):
		rdf_node = self.nspace[nodeId]
		#print rdf_node
		node_triples = []
		if attr_type == "set" or attr_type == "list":
			for term in nodeInfo[attr]:
				pr_term = condition(term) if condition else term
				node_triples.append((rdf_node, pred, pr_term))
		else:
			pr_term = condition(nodeInfo[attr]) if condition else nodeInfo[attr]
			print pr_term
			node_triples.append((rdf_node, pred, pr_term))
		return node_triples

	def gen_edge_triples(self, edge, edgeInfo):
		#print edge
		edge_triples = []
		edge_node = BNode()
		edge_type_f = lambda x: self.nspace[x.title()]
		pr_type = edge_type_f(edgeInfo["type"]) if "type" in edgeInfo else self.nspace.Resource
		edge_triples.append((edge_node, RDF.type, pr_type)) # assumes each edge is typed using "type"
		rdf_node1 = self.nspace[edge[0]]
		rdf_node2 = self.nspace[edge[1]]
		edge_triples.append((edge_node, self.nspace.source, rdf_node1))
		edge_triples.append((edge_node, self.nspace.target, rdf_node2))
		return edge_triples

	def save_rdfG(self, filename):
		rdf_file = open(filename, "w+")
		for l in self.rdfG.serialize(format='nt').splitlines():
			print l
			rdf_file.write(l.decode('ascii') + "\n")
		rdf_file.close()
		return True

def testing_classes():
	'''This includes simple testing functions'''
	G = sys.argv[1]
	#presenter = GraphPresenter(G)
	#rdfG = presenter.gen_rdfG("http://onto-apps.stanford.edu/phlegra/")
	#presenter.save_rdfG("drug-network.nt")
	manipulator = GraphManipulator(G)
	drug_node = "drug_195"
	dise_node = "dise_9758"
	

if __name__ == "__main__":
	testing_classes()
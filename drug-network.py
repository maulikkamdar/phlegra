import hydra.tpf
import rdflib
import networkx as nx
import json, re, argparse
from collections import Counter
from copy import deepcopy
#import flask

# This code is a mess ....
# @TODO umm, a lot of changes, include argparse, clean up, deal seamlessly with auxiliary functions to post-process data, etc. - Check for errors in Bio2RDF ...
# This code is the gist of everything ... It queries DrugBank, KEGG, Comparative Toxicogenomics Database, and PharmGKB
# The RDF datasets are first converted to HDT and exposed through Triple Pattern Fragment Servers
# Ideally this should serve as a federated query system, where we are given a list of SPARQL queries and configuration details
# The problem is that there are several issues in the Bio2RDF converted RDF stores, that hence require all the other auxiliary functions
# The ideal query would be a CONSTRUCT query that retrieves all the SELECT data from the different sources, and generates a Drug-Protein-Pathway-Phenotype network, 
# that can be stored as a RDF graph available for query.
# However, this code currently creates a NetworkX graph (that needs to be again translated to an RDF graph .... an extra step)
# Because of the auxiliary post-processing of the SELECT query data and SPARQL endpoint limitations ....

datasources = ["drugbank", "kegg", "cas", "atc", "pharmgkb", "sider", "mesh", "ctd", "hgnc", "genbank", "gi", "uniprot", "ncbigene", "omim", "icd10", "go"]
# this is hardcoded for a single term
#conterm = "<http://bio2rdf.org/drugbank:DB00203>"

config_details = read_json("config-default.json")
prefix_header = ""
for prefix in config_details["prefixes"]:
	prefix_header = prefix_header + "PREFIX " + prefix + ": <" + config_details["prefixes"][prefix] + "> " 

queries = read_json("sparql-queries.json")

G = nx.DiGraph()
ids = {"drugs": {"drugbank_ids": {}, "cas_ids": {}, "kegg_ids": {}, "pharmgkb_ids": {}, "dailymed_ids": {}, "identifers": {}, "mesh_ids": {}}, 
"genes": {"hgnc_ids": {}, "drugbank_ids": {}, "gi_ids": {}, "genbank_ids":{}, "uniprot_ids":{}, "identifers": {}, "kegg_ids": {}, "pharmgkb_ids": {}, "ncbigene_ids": {}, "omim_ids": {}}, 
"dises": {"omim_ids": {}, "mesh_ids": {}, "meddra_ids":{}, "umls_ids":{}, "identifers": {}, "icd10_ids": {}, "kegg_ids": {}}, "paths": {"kegg_ids": {}, "identifers": {}}, "goterms": {"identifers": {}, "go_ids": {}}}


def read_json(filename):
	with open(filename) as config:
		details = json.load(config)
	return details

def getPredicate(term):
	predicateParts = re.split('[^0-9a-zA-Z]+', term[1:-1])
	return predicateParts[len(predicateParts)-1]

def getId(uri):
	uriParts = uri[1:-1].split(":")
	return uriParts[len(uriParts)-1]

def get_query(query, term=None, datasource="drugbank"):
	query = prefix_header + "".join(queries[datasource][query])
	if term:
		query = query.replace("TERM_URI", term)
	print query
	return query

def get_endpoint(hdtname):
	g = rdflib.Graph('TPFStore')
	g.open(config_details["endpoints"][hdtname]["link"])
	return g

def get_spo(row):
	subj = row[0].n3()
	pred = getPredicate(row[1].n3())
	obj = row[2].n3()
	return (subj, pred, obj)

def extract_list(datasource, hdtname, ordering, datatype="drugs", altdb=None):
	g = get_endpoint(hdtname)
	query = get_query("all"+datatype, datasource=hdtname)
	results = g.query(query)
	for row in results:
		#because of URI issue in Pharmgkb x-refs
		if hdtname == "pgkbgenes":
			row = list(row)
			for k in range(len(row)):
				if ordering[k] != "pharmgkb" and ordering[k] != "title" and row[k]:
					row[k] = rdflib.term.URIRef(row[k][:-1])
		term_id = getId(row[0].n3())
		#because of orthologous genes in kegg
		if datasource == "kegg" and term_id[0:2] == "ko":
			continue
		if not term_id in ids[datatype][datasource + "_ids"]:
			identifer = datatype[:-1] + "_" + str(len(ids[datatype]["identifers"]))
			source = altdb if altdb else datasource
			G.add_node(identifer, title=set([]), type=datatype, source=set([source]))
			ids[datatype][datasource + "_ids"][term_id] = identifer
			ids[datatype]["identifers"][identifer] = {}
			for source in datasources:
				ids[datatype]["identifers"][identifer][source] = set([])
			ids[datatype]["identifers"][identifer][ordering[0]].add(term_id)
		else:
			identifer = ids[datatype][datasource + "_ids"][term_id]
			source = altdb if altdb else datasource
			G.node[identifer]["source"].add(source)
		if row[1]:
			G.node[identifer]["title"].add(row[1].n3())
		for k in range(2, len(row)):
			if row[k]:
				if len(row[k].n3()) == 0:
					continue
				db_id = getId(row[k].n3())
				if ordering[k] == "kegg":
					db_id = db_id.upper()
				if db_id not in ids[datatype][ordering[k] + "_ids"]:
					ids[datatype][ordering[k] + "_ids"][db_id] = identifer
					ids[datatype]["identifers"][identifer][ordering[k]].add(db_id)
	print "Extracted " + datatype + " List from " + datasource


def add_net_edge(resource1, resource2, datatype1, datatype2, datasource, pred, assoc=None, datasource2=None):
	#print (resource1, resource2, datatype1, datatype2, datasource, pred, assoc, datasource2)
	if datasource2:
		if not resource1 in ids[datatype1][datasource + "_ids"] or not resource2 in ids[datatype2][datasource2 + "_ids"]:
			return None
		else:
			id2 = ids[datatype2][datasource2 + "_ids"][resource2]
	else:
		if not resource1 in ids[datatype1][datasource + "_ids"] or not resource2 in ids[datatype2][datasource + "_ids"]:
			return None
		else:
			id2 = ids[datatype2][datasource + "_ids"][resource2]
	id1 = ids[datatype1][datasource + "_ids"][resource1]
	
	if G.has_edge(id1, id2):
		G[id1][id2]["source"].append(datasource)
		if assoc:
			G[id1][id2]["assoc"].append(assoc)
	else:
		if assoc:
			G.add_edge(id1, id2, type=pred, source=[datasource], assoc=[assoc])
		else:
			G.add_edge(id1, id2, type=pred, source=[datasource], assoc=[])

def modify_resource_2(datasource, term):
	resources = []
	if datasource == "kegg" and len(term.split()) > 1:
		objParts = term.split()
		namespaceParts = objParts[0].split("_")
		if len(namespaceParts) > 1:
			for k in range(len(objParts)):
				if k == 0:
					resource = getId(objParts[k])
					resources.append(resource)
				else:
					resources.append(getId(namespaceParts[0] + objParts[k]))
	elif datasource == "ncbigene":
		resource = "map"+ getId(term.n3())
		resources.append(resource)
	else:
		resource = getId(term.n3())
		#if resource[0] != "K":
		resources.append(resource)
	return resources

def extract_edges(datasource, hdtname, query_name, pred, datasource2=None):
	g = get_endpoint(hdtname)
	query = get_query(query_name, datasource=hdtname)
	results = g.query(query)
	for row in results:
		resource_1 = getId(row[0].n3())
		if datasource == "kegg" and pred == "pathway":
			resource_1 = resource_1.upper()
		elif datasource == "kegg" and pred == "disease":
			resource_1 = "map" + resource_1
		resources = modify_resource_2(datasource, row[1])
		for resource_2 in resources:
			assoc = []
			for k in range(2, len(row)):
				assoc.append(row[k].n3())
			if pred == "interacts":
				add_net_edge(resource_1, resource_2, "drugs", "drugs", datasource, pred, assoc, datasource2=datasource2)
			elif pred == "target":
				add_net_edge(resource_1, resource_2, "drugs", "genes", datasource, pred, assoc, datasource2=datasource2)
			elif pred == "pathway":
				add_net_edge(resource_1, resource_2, "genes", "paths", datasource, pred, assoc, datasource2=datasource2)
			elif pred == "disease":
				add_net_edge(resource_1, resource_2, "paths", "dises", datasource, pred, assoc, datasource2=datasource2)
			elif pred in ["component", "function", "process"]:
				#print (resource_1, resource_2, "genes", "goterms", "uniprot", pred)
				add_net_edge(resource_1, resource_2, "genes", "goterms", "uniprot", pred, assoc, datasource2=datasource2)
			else:
				add_net_edge(resource_2, resource_1, "genes", "drugs", datasource, pred, assoc)

def extract_edges_red(datasource, hdtname, query_name, pred, datasource2=None):
	g = get_endpoint(hdtname)
	query = get_query(query_name, datasource=hdtname)
	results = g.query(query)
	for row in results:
		#print row
		resource_1 = getId(row[0].n3())
		entities = resource_1.split("-")
		id1 = ids["drugs"]["mesh_ids"][entities[0]]
		id2 = ids["genes"]["ncbigene_ids"][entities[1]]
		if G.has_edge(id1, id2):
			G[id1][id2]["assoc"].append(row[1].n3())
			G[id1][id2]["source"].append("ctd")
		else:
			G.add_edge(id1, id2, type=pred, source=["ctd"], assoc=[row[1].n3()])


def output_ids(filename):
	print_ids = deepcopy(ids)
	for datatype in ids:
		for identifer in ids[datatype]["identifers"]:
			for source in ids[datatype]["identifers"][identifer]:
				print_ids[datatype]["identifers"][identifer][source] = list(ids[datatype]["identifers"][identifer][source])
	with open(filename, "w+") as f:
		json.dump(print_ids, f)
	f.close()

def read_ids(filename):
	ids = read_json(filename)
	for datatype in ids:
		for identifer in ids[datatype]["identifers"]:
			for source in ids[datatype]["identifers"][identifer]:
				ids[datatype]["identifers"][identifer][source] = set(ids[datatype]["identifers"][identifer][source])
	return ids	

def print_graph(filename):
	print_G = deepcopy(G)
	for node in G.nodes():
		if G.node[node]['type'] == "goterms":
			print_G.remove_node(node)
		else:
			sources = list(G.node[node]["source"])
			sources = [x for x in sources if x is not None]
			print_G.node[node]["title"] = "; ".join(list(G.node[node]["title"]))
			print_G.node[node]["source"] = "; ".join(sources)
	for edge in G.edges():
		print G[edge[0]][edge[1]]
		print_G[edge[0]][edge[1]]["source"] = "; ".join(list(G[edge[0]][edge[1]]["source"]))
		print_G[edge[0]][edge[1]]["assoc"] = "; ".join(["-".join(k) for k in G[edge[0]][edge[1]]["assoc"]])
	nx.write_graphml(print_G, filename)


def get_gene_int_id(term, actnot):
	if term in ids["genes"][actnot+"_ids"]:
		identifer = ids["genes"][actnot+"_ids"][term]
		return identifer
	else:
		return None

def get_gene_int(iterpart):
	gene_node = None
	altdbIds = iterpart.split("|")
	for altdbid in altdbIds:
		altid = altdbid.split(":")
		if altid[0] == "uniprotkb":
			gene_node = get_gene_int_id(altid[1], "uniprot")
			if gene_node:
				break
		elif altid[0] == "genbank_protein_gi":
			gene_node = get_gene_int_id(altid[1], "gi")
			if gene_node:
				break
			gene_node = get_gene_int_id(altid[1], "genbank")
			if gene_node:
				break
		elif altid[0] == "hgnc":
			gene_node = get_gene_int_id(altid[1], "hgnc")
			if gene_node:
				break
	return gene_node

def overlay_interactions(filename):
	irefindex = open(filename)
	interactions = irefindex.readlines()
	irefindex.close()
	print len(interactions)
	rolesFile = open("roles.tsv")
	roles = rolesFile.readlines()
	rolesFile.close()
	print len(roles)
	rolesMatrix = {}
	for k in range(len(roles)):
		if k == 0:
			continue
		rolesParts = roles[k].strip().split("\t")
		if not rolesParts[0] in rolesMatrix:
			rolesMatrix[rolesParts[0]] = {}
		rolesMatrix[rolesParts[0]][rolesParts[1]] = {"assoc": rolesParts[3], "source": rolesParts[2]}
		
	def add_gene_int_edge(id1, id2, pred, datasource, assoc):
		if G.has_edge(id1, id2):
			G[id1][id2]["source"].append(datasource)
			G[id1][id2]["assoc"].append(assoc)
		else:
			G.add_edge(id1, id2, type=pred, source=[datasource], assoc=[assoc])
			
	for k in range(len(interactions)):
		iterParts = interactions[k].strip().split("\t")
		print len(iterParts)
		gene_node1 = None
		gene_node2 = None
		if not gene_node1:
			gene_node1 = get_gene_int(iterParts[0])
		if not gene_node2:
			gene_node2 = get_gene_int(iterParts[1])
		if not gene_node1:
			gene_node1 = get_gene_int(iterParts[2])
		if not gene_node2:
			gene_node2 = get_gene_int(iterParts[3])
		if not gene_node1:
			gene_node1 = get_gene_int(iterParts[4])
		if not gene_node2:
			gene_node2 = get_gene_int(iterParts[5])
		if gene_node1 and gene_node2:
			print (gene_node1, gene_node2)
			role1 = iterParts[8]
			role2 = iterParts[9]
			source_node = rolesMatrix[role1][role2]["source"]
			if source_node == "gene1":
				add_gene_int_edge(gene_node1, gene_node2, "interacts", iterParts[7], rolesMatrix[role1][role2]["assoc"])
			else:
				add_gene_int_edge(gene_node2, gene_node1, "interacts", iterParts[7], rolesMatrix[role1][role2]["assoc"])

def output_names(filename, datatype, idtype=None):
	f = open(filename, "w+")
	for node in G.nodes():
		if G.node[node]["type"] == datatype:
			title = "; ".join(list(G.node[node]["title"]))
			title = title.encode('utf-8').strip()
			print node, title
			if idtype:
				identifer = ""
				if idtype in ids[datatype]["identifers"][node]:
					identifer = ";".join(list(ids[datatype]["identifers"][node][idtype]))
				f.write("%s\t%s\t%s\n" % (node, title, identifer))
			else:
				f.write("%s\t%s\n" % (node, title))
	f.close()

def ctd_genes_dises(filename):
	gdisindex = open(filename)
	interactions = gdisindex.readlines()
	gdisindex.close()
	print len(interactions)
	for k in range(len(interactions)):
		inter = interactions[k].strip().split("\t")
		gene_node = None
		dise_node = None
		mesh_parts = inter[2].split(":")
		if len(mesh_parts) == 1:
			continue
		mesh_term = mesh_parts[1]
		if inter[0] in ids["genes"]["ncbigene_ids"]:
			gene_node = ids["genes"]["ncbigene_ids"][inter[0]]
		if mesh_term in ids["dises"]["mesh_ids"]:
			dise_node = ids["dises"]["mesh_ids"][mesh_term]
		if gene_node and dise_node:
			G.add_edge(gene_node, dise_node, type="disease", source=["ctd"], assoc=[])

def main():
	extract_list("drugbank", "drugbank", ("drugbank", "title", "kegg", "cas", "pharmgkb"), "drugs")
	extract_list("drugbank", "drugbank", ("drugbank", "title", "hgnc", "uniprot", "genbank", "gi"), "genes")
	extract_edges("drugbank", "drugbank", "drugtargets", "target")
	extract_edges("drugbank", "drugbank", "drugenzymes", "enzyme")
	extract_edges("drugbank", "drugbank", "drugtransporters", "transporter")

	extract_list("kegg", "kdrugs", ("kegg", "title", "drugbank", "cas", "pharmgkb"), "drugs")
	extract_list("hgnc", "kgenes", ("hgnc", "title", "gi", "uniprot", "omim", "kegg", "ncbigene"), "genes", "kegg")
	extract_edges("kegg", "kdrugs", "drugtargets", "target")
	extract_edges("kegg", "kdrugs", "drugenzymes", "enzyme")

	extract_list("pharmgkb", "pgkbdrugs", ("pharmgkb", "title", "drugbank", "cas", "kegg"), "drugs")
	extract_list("kegg", "kpaths", ("kegg", "title"), "paths")
	extract_edges("kegg", "kgenes", "genepaths", "pathway")
	extract_list("kegg", "kdises", ("kegg", "title", "icd10", "mesh", "omim"), "dises")
	extract_edges("kegg", "kpaths", "pathdises", "disease")

	extract_list("hgnc", "pgkbgenes", ("hgnc", "title", "uniprot", "pharmgkb", "omim", "ctd"), "genes", "pharmgkb")
	extract_edges("pharmgkb", "pgkbrel", "druggenes", "target")

	extract_list("go", "go", ("go", "title"), "goterms")
	extract_edges("uniprot", "goa", "component", "component", datasource2="go")
	extract_edges("uniprot", "goa", "function", "function", datasource2="go")
	extract_edges("uniprot", "goa", "process", "process", datasource2="go")

	extract_list("cas", "ctddrugs", ("cas", "title", "mesh"), "drugs", "ctd")
	extract_list("pharmgkb", "ctdgenes", ("pharmgkb", "title", "ncbigene", "uniprot"), "genes", "ctd")
	extract_list("mesh", "ctddises", ("mesh", "title"), "dises", "ctd")

	extract_edges_red("mesh", "redctddruggenes", "druggenes", "target", datasource2="ncbigene") # Query the reduced CTD data source, as the CTD data source exposed through TPF server time outs due to size of CTD - not sure if this is a problem of Hydra or TPF :(
	extract_edges("ncbigene", "ctdgenepaths", "genepaths", "pathway", datasource2="kegg")
	extract_edges("kegg", "ctdpathdises", "pathdises", "disease", datasource2="mesh")

	#overlay_interactions("../data/irefindex.tsv")
	G.remove_nodes_from(nx.isolates(G))
	nx.write_gpickle(G, "drug-network.gpickle")
	output_names("drug-names.tsv", "drugs", "mesh")

if __name__ == "__main__":
	main()



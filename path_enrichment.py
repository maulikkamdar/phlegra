import networkx as nx
from utils import MatrixIO
from copy import deepcopy

G = nx.read_gpickle("reducedGraph.gpickle")
ones_pairs_file = open("baseline_pair_1s.tsv")

def gen_keys():
	key_file = open("fmap_key_file.json")
	key_lines = key_file.readlines()
	key_file.close()
	keys = {}
	for k in key_lines:
		key_parts = k.strip().split("\t")
		key_type = key_parts[0][0:4]
		if key_type == "drug":
			key_id = "dr" + key_parts[1]
		elif key_type == "dise":
			key_id = "di" + key_parts[1]
		elif key_type == "gene":
			key_id = "ge" + key_parts[1]
		elif key_type == "path":
			key_id = "pa" + key_parts[1] 
		keys[key_parts[0]] = key_id
	return keys

#keys = gen_keys()
ones_pairs = ones_pairs_file.readlines()
ones_pairs_file.close()
print "read graph"

def create_pathlib():
	path_lib = {}
	for edge in G.edges():
		path_lib[edge[0] + "->" + edge[1]] = 0
		path_lib[edge[1] + "->" + edge[0]] = 0
	return path_lib

path_lib = create_pathlib()

gene_paths = [n for n in G.nodes() if G.node[n]["type"] in ["genes", "paths"]]
assoc_paths = {}
no_paths = []

mfio = MatrixIO()
pc_count = open("pc_count.tsv", "w+")

for k in range(len(ones_pairs)):
	pr_parts = ones_pairs[k].strip().split()
	dreg = deepcopy(gene_paths)
	dreg.extend([pr_parts[0], pr_parts[1]])
	all_paths = []
	assoc_id = pr_parts[0] + "-->" + pr_parts[1]
	subG = G.subgraph(dreg)
	rel_nodes = nx.shortest_path(subG, pr_parts[0]).keys()
	if k%1000 == 0:
		mfio.save_matrix(path_lib, "path_lib.dat")
		mfio.save_matrix(assoc_paths, "assoc_paths_" + str(k) + ".dat")
		mfio.save_matrix(no_paths, "no_paths_" + str(k) + ".dat")
		assoc_paths = {}
		no_paths = []
	if not pr_parts[1] in rel_nodes:
		no_paths.append(assoc_id)
		continue
	for m in nx.all_shortest_paths(subG, pr_parts[0], pr_parts[1]):
		for n in range(len(m)-1):
			path = m[n] + "->" + m[n+1]
			path_lib[path] += int(pr_parts[2])
		#all_paths.append("->".join([keys[mpart] for mpart in m]))
		all_paths.append("->".join(m))
	assoc_paths[assoc_id] = all_paths
	pc_count.write(str(k) + "\t" + str(len(all_paths)) + "\n")
	#print str(k) + "\t" + str(len(all_paths))

pc_count.close()

mfio.save_matrix(path_lib, "path_lib.dat")
mfio.save_matrix(assoc_paths, "assoc_paths_" + str(k)+ ".dat")
mfio.save_matrix(no_paths, "no_paths_" + str(k) + ".dat")
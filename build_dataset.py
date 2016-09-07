import os, sys, argparse
import numpy as np
import networkx as nx

# This script rewrites the FAERS data files generated from parseFAERS.py to standard ids used in the drug-protein-path-disease network
# Input: Network generated by drug-network.py, Folder of parsed FAERS files, mapped_locator files (generated from reaction-mapper.py) and mapped drug files (generated from drug-term.py)
# Output: Key file, that maps ids to different ids used for different entities (e.g. drugbank, pharmgkb ids for drugs)
# Output: Transformed data file that maps patient id to drug ids
# Output: Transformed labels file that maps patient id to reaction ids
# PS. I have no idea what the below two lines of code did previously --- was I using the strongest connected component and if yes, why did I stop using it -- did I save the strongest connected component separately under drug-network.py (That line does not exist somehow!)
# python build_dataset.py -g to_use_data_files/drug_network_reduced2.gpickle -faers /Users/mkamdar/Desktop/PhD/data/fda/fda/ -drugpref to_use_data_files/map_drug -reactpref to_use_data_files/nmap -indipref to_use_data_files/nmap -o to_use_data_files/fmap

#Gc = nx.Graph(nx.read_gpickle("drug-gene-path.gpickle")) # read_model
#G = max(nx.connected_component_subgraphs(Gc), key=len) # work on the larges connected component only 

keys = {"drugs": {}, "dises": {}, "genes": {}, "paths": {}}
kount = {"drugs":0, "dises":0, "genes":0, "paths":0}

def parse_args():
	p = argparse.ArgumentParser(description='Builds the Data file and the Reaction/Indications (Labels files)')
	p.add_argument('-g', '--graph', help='Location of the drug-protein-path-disease network generated by drug-network.py', default="drug_network_reduced.gpickle")
	p.add_argument('-faers', '--faers', help='Folder of the FAERS file generated by parseFAERS.py', default="../data/fda")
	p.add_argument('-drugpref', '--drugpref', help='Prefix of the Drug mapping files', default="results")
	p.add_argument('-reactpref', '--reactpref', help='Prefix of the Reaction mapping files', default="results")
	p.add_argument('-indipref', '--indipref', help='Prefix of the Indication mapping files', default="results")
	p.add_argument('-o', '--prefix', help='prefix of the output file (with directory)', default="results")
	return p.parse_args()

def generate_drug_dict(drug_pref):
	drugs = {}
	for k in range(1, 3):
		drug_file = open(drug_pref + "_match_" + str(k) + ".tsv")
		drug_lines = drug_file.readlines()
		drug_file.close()
		for m in range(len(drug_lines)):
			drug_parts = drug_lines[m].strip().split("\t")
			main_drug = drug_parts[2]
			if main_drug in keys["drugs"]:
				drugs[drug_parts[1].strip()] = [main_drug]
				print main_drug
			else:
				main_drug = None
			if k == 2 and drug_parts[3] in keys["drugs"]:
				if main_drug:
					drugs[drug_parts[1].strip()].append(drug_parts[3])
				else:
					drugs[drug_parts[1].strip()] = [drug_parts[3]]
	return drugs

def generate_phenotype_dict(phen_pref, phen_type):
	#phen_pref + "_" + phen_type + "_locatorfile1.tsv", "w+"
	phenotypes = {}
	for k in range(0,2):
		reaction_file = open(phen_pref + "_" + phen_type + "_locatorfile" + str(k) + ".tsv")
		reactionlines = reaction_file.readlines()
		reaction_file.close()
		for m in range(len(reactionlines)):
			reaction_parts = reactionlines[m].strip().split("\t")
			if k == 0:
				reaction = reaction_parts[2]
			else:
				reaction = reaction_parts[1]
			if reaction in keys["dises"]:
				phenotypes[reaction_parts[0].strip()] = reaction
	return phenotypes

def populate_keys(key_file_name):
	pass

def main():
	args = parse_args()
	G = nx.read_gpickle(args.graph)
	key_file_name =  args.prefix + "_key_file.json" 
	data_file_name = args.prefix + "_data_file.tsv"
	reaction_file_name = args.prefix + "_reaction_labels.tsv"
	indication_file_name = args.prefix + "_indication_labels.tsv"

	print "working with model of " + str(len(G.nodes())) + " nodes and " + str(len(G.edges())) + " edges"

	keyfile = open(key_file_name, "w+")

	for node in sorted(G.nodes()):
		type = G.node[node]["type"]
		keys[type][node] = kount[type]
		keyfile.write(node + "\t" + str(kount[type]) + "\n")
		kount[type] = kount[type] + 1

	keyfile.close()
	
	drugs = generate_drug_dict(args.drugpref)
	reactions = generate_phenotype_dict(args.reactpref, "reaction")
	indications  = generate_phenotype_dict(args.indipref, "indication")
	data = {}

	print "output file"
	print "Total number of drugs:"
	print len(drugs)

	print "Total number of reactions:"
	print len(reactions)

	print "Total number of indications:"
	print len(indications)

	data_file = open(data_file_name, "w+")
	reactions_file = open(reaction_file_name, "w+")
	indications_file = open(indication_file_name, "w+")
	folder_path = args.faers
	num_patients = 0
	
	for root, dirs, files in os.walk(folder_path):
		for file in files:
			filename = file.split(".")
			if filename[1] != "tsv":
				continue
			print file
			aers_file = open(folder_path + file)
			aers_lines = aers_file.readlines()
			aers_file.close()
			for k in range(len(aers_lines)):
				aers_parts = aers_lines[k].strip().split("\t")
				aers_drugs = aers_parts[1].split(";")
				aers_reactions = aers_parts[2].split(";")
				if len(aers_parts) > 3:
					aers_indications = aers_parts[3].split(";")
				else:
					aers_indications = []
				model_drugs = []
				model_reactions = []
				model_indications = []

				for drug in aers_drugs:
					if drug.strip() in drugs:
						model_drugs.extend(drugs[drug.strip()])
				for reaction in aers_reactions:
					if reaction.strip() in reactions:
						model_reactions.append(reactions[reaction.strip()])
				for indication in aers_indications:
					if indication.strip() in indications:
						model_indications.append(indications[indication.strip()])

				if len(model_drugs) > 0 and len(model_reactions) > 0:
					num_patients = num_patients + 1
					print model_drugs
					print model_reactions
					print model_indications
					for drug in sorted([keys["drugs"][drug] for drug in model_drugs]):
						data_file.write(str(aers_parts[0]) + "\t" + str(drug) + "\n")
					for dise in sorted([keys["dises"][dise] for dise in model_reactions]):
						reactions_file.write(str(aers_parts[0]) + "\t" + str(dise) + "\n")
					for dise in sorted([keys["dises"][dise] for dise in model_indications]):
						indications_file.write(str(aers_parts[0]) + "\t" + str(dise) + "\n")

	print num_patients
	data_file.close()
	reactions_file.close()
	indications_file.close()

if __name__ == "__main__":
	main()


















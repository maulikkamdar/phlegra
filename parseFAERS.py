from lxml import etree
from xml.sax.saxutils import unescape
import hashlib, os, sys

# This script parses the XML files of FAERS to extract only the drug names, adverse reactions and drug indications and convert the XML files to TSV files with 4 columns - Universal patient No (auto index), drug names, adverse reactions and drug indications.
# This script also publishes the set of unique drug names and unique reaction/indication names that can be later used to associate RxNorm or MEDDRA terms 
# Input: folder_name where the FDA files are stored, prefix for the named files of drugs and reactions, indications
# Output: Named file for Drugs, Reactions, and Indications as well as Tsv files are stored in the same folder
# python parseFAERS.py FAERS_folder DRUG_file REACTION_file INDICATION_file
# e.g. python parseFAERS.py /Users/mkamdar/Desktop/PhD/data/fda/ "drug-faers.tsv" "reaction-faers.tsv" "indications-faers.tsv"
# @TODO use argparse instead of sys
# @TODO to include a method to store reactions and indications separately so that we can map them to different UMLS terminologies (MEDDRA, SNOMED CT)
# @TODO try to use standardized FAERS mining scripts published in the 2016 Banda et al. paper

def get_reaction(node):
	reaction = None
	for tag in range(len(node)):
		if node[tag].tag == "reactionmeddrapt":
			reaction = node[tag].text
	return reaction

def get_drug(node):
	''' This function also returns indications, as FAERS reports stores indication and drug name in the same node
	'''
	indication = None
	drug = None
	for tag in range(len(node)):
		if node[tag].tag == "drugindication":
			indication = node[tag].text # is it a chance that the same drug can be given for multiple indications?
			#indications.add(indication)
		elif node[tag].tag == "medicinalproduct":
			drug = node[tag].text
		elif node[tag].tag == "activesubstance":
			for m in range(len(node[tag])):
				if node[tag][m].tag == "activesubstancename":
					drug = node[tag][m].text
	return (indication, drug)

def print_file(file_info, entity_info):
	entities = sorted(list(entity_info))
	for entity in entities:
		file_info.write(entity + "\n")
	file_info.close()

def main():
	drugs = set([])
	reactions = set([])
	indications = set([])
	patient_no = 0
	folder_path = sys.argv[1]
	prefix = sys.argv[2]
	drugfile = open(prefix + "_drugs.tsv", "w+")
	reactionfile = open(prefix + "_reactions.tsv", "w+")
	indicationfile = open(prefix + "_indications.tsv", "w+")

	for root, dirs, files in os.walk(folder_path):
		for file in files:
			filename = file.split(".")
			if filename[1] != "xml":
				continue
			print file
			data_file = open(folder_path + "faers_" + filename[0] + ".tsv", "w+")
			xmldoc = etree.parse(folder_path + file)
			for elt in xmldoc.getiterator("safetyreport"):
				for patient in elt.getiterator("patient"):
					patient_info = {"reactions":set([]), "drugs": set([]), "indication": set([])}
					for tag in range(len(patient)):
						if patient[tag].tag == "reaction":
							reaction_info = get_reaction(patient[tag])
							if reaction_info:
								reactions.add(reaction_info)
								patient_info["reactions"].add(reaction_info)
						elif patient[tag].tag == "drug":
							drug_info = get_drug(patient[tag])
							if drug_info[1]:
								drugs.add(drug_info[1])
								patient_info["drugs"].add(drug_info[1])
							if drug_info[0]:
								indications.add(drug_info[0])
								patient_info["indication"].add(drug_info[0])

					data_file.write("%s\t%s\t%s\t%s\n" % (patient_no, "; ".join(list(patient_info["drugs"])), "; ".join(list(patient_info["reactions"])), "; ".join(list(patient_info["indication"]))))
					patient_no +=1
			data_file.close()

	print_file(drugfile, drugs)
	print_file(reactionfile, reactions)
	print_file(indicationfile, indications)

if __name__ == '__main__':
	main()



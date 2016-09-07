import os, sys, argparse
import operator

# This script generates 2 files - the first file has 3 columns: Drug, Adverse Event, Count
# The second file has 3 columns: Drug-Drug Interaction, Adverse Event, Count
# Input: Drugs File, Reactions file generated from build_dataset.py
# Output: the 2 files, that will be fed into the R script to generate the different measures
# python generate_counts.py -keyfile to_use_data_files/fmap_key_file.json -drugfile to_use_data_files/fmap_data_file.tsv -reactfile to_use_data_files/fmap_reaction_labels.tsv -o to_use_data_files/baseline

pair_1s = {}
pair_2s = {}

def parse_args():
	p = argparse.ArgumentParser(description='Generate count statistics for 1-1 associations and 2-1 associations')
	p.add_argument('-keyfile', '--keyfile', help='Location of the drug-gene-path-dise key file generated from build_dataset.py', default="key_file.json") # This is actually a TSV file
	p.add_argument('-drugfile', '--drugfile', help='Location of the data_file', default="data_file.tsv")
	p.add_argument('-reactfile', '--reactfile', help='Location of the reactions_file', default="reactions_file.tsv")
	#p.add_argument('-N', '--N', help="total number of patients", default="3037542") # this should be integer
	p.add_argument('-o', '--prefix', help='prefix of the output file (with directory)', default="results")
	return p.parse_args()

def populate_patient_dict(file_name):
	patient_dict = {}
	patient_file = open(file_name)
	patient_lines = patient_file.readlines()
	patient_file.close()

	for k in range(len(patient_lines)):
		patient_parts = patient_lines[k].strip().split()
		if not patient_parts[0] in patient_dict:
			patient_dict[patient_parts[0]] = []
		patient_dict[patient_parts[0]].append(patient_parts[1])

	return patient_dict

def gen_1_pairs(drugs, dises):
	pairs = set([])
	for drug in drugs:
		for dise in dises:
			pairs.add((drug, dise))
	return pairs

def gen_2_pairs(drugs, dises):
	pairs = set([])
	for drug1 in drugs:
		for drug2 in drugs:
			if int(drug1) < int(drug2):
				for dise in dises:
					pairs.add(((drug1, drug2), dise))
	return pairs

def main():
	args = parse_args()
	key_file = open(args.keyfile)
	dise_dict = {}
	drug_dict = {}

	for k in key_file.readlines():
		key_params = k.strip().split()
		if key_params[0][0:4] == "dise":
			dise_dict[key_params[1]] = key_params[0]
		elif key_params[0][0:4] == "drug":
			drug_dict[key_params[1]] = key_params[0]
	key_file.close()

	print "Generated " + str(len(drug_dict)) + " Drugs"
	print "Generated " + str(len(dise_dict)) + " Diseases"

	drug_list = populate_patient_dict(args.drugfile)
	react_list = populate_patient_dict(args.reactfile)
	print "Saved Drug infor on "+ str(len(drug_list)) + " Patients"
	print "Saved React infor on "+ str(len(react_list)) + " Patients"

	for patient in drug_list:
		if patient in react_list:
			pair_1d = gen_1_pairs(drug_list[patient], react_list[patient])
			for pair in pair_1d:
				if not pair in pair_1s:
					pair_1s[pair] = 0
				pair_1s[pair] += 1
			if len(drug_list[patient]) > 1:
				pair_2d = gen_2_pairs(drug_list[patient], react_list[patient])
				for pair in pair_2d:
					if not pair in pair_2s:
						pair_2s[pair] = 0
					pair_2s[pair] += 1

	sorted_pair1s = sorted(pair_1s.items(), key=operator.itemgetter(1), reverse=True)
	sorted_pair2s = sorted(pair_2s.items(), key=operator.itemgetter(1), reverse=True)
	
	def save_dict(pair_dict, file_name, dual=False):
		dict_file = open(file_name, "w+")
		for pair in pair_dict:
			print pair
			if dual:
				entry = (drug_dict[pair[0][0][0]], drug_dict[pair[0][0][1]])
			else:
				entry = drug_dict[pair[0][0]]
			dict_file.write(str(entry) + "\t" + dise_dict[pair[0][1]] + "\t" + str(pair[1]) + "\n")
		dict_file.close()

	save_dict(sorted_pair1s, args.prefix + "_pair_1s.tsv")
	save_dict(sorted_pair2s, args.prefix + "_pair_2s.tsv", dual=True)

if __name__ == "__main__":
	main()
	

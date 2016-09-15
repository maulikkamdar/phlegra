import os, sys, argparse
import operator

# This script generates 2 files - the first file has expected counts for 1s pairs and second file has expected counts for the 2s pairs
# It also generates auxiliary files on drug count information
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
	p.add_argument('-pair1s', '--pair1s', help='Location of the baseline 1s pairs', default="baseline_pairs_1s.tsv")
	# add here another argument to include 2s pairs
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

def generate_counts(patient_dict, id_dict, count_filename):
	count_dict = {}
	for patient in patient_dict:
		for k in patient_dict[patient]:
			k_id = id_dict[k]
			if not k_id in count_dict:
				count_dict[k_id] = 0
			count_dict[k_id] += 1
	
	sorted_x = sorted(count_dict.items(), reverse=True, key=operator.itemgetter(1))
	count_file = open(count_filename, "w+")
	for k in range(len(sorted_x)):
		count_file.write(sorted_x[k][0] + "\t" + str(sorted_x[k][1]) + "\n")
	count_file.close()
	return count_dict


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

	drug_count = generate_counts(drug_list, drug_dict, "drug_count.tsv")
	react_count = generate_counts(react_list, dise_dict, "dise_count.tsv")

	print "Generated " + str(len(drug_count)) + " drugs"
	print "Generated " + str(len(react_count)) + " reactions"

	pair_file = open(args.pair1s)
	pair_lines = pair_file.readlines()
	pair_file.close()

	exp_file = open(args.prefix + "_pair1s.tsv", "w+")

	for k in range(len(pair_lines)):
		pair_params = pair_lines[k].strip().split("\t")
		# move this to deal with dual interactions
		if pair_params[0] in drug_count and pair_params[1] in react_count:
			expcount = float(drug_count[pair_params[0]]*react_count[pair_params[1]])/len(drug_list)
		exp_file.write(pair_params[0] + "\t" + pair_params[1] + "\t" + str(expcount) + "\n")

	exp_file.close()

if __name__ == "__main__":
	main()	

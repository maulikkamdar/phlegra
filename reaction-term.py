import os, json, argparse

# FAERS reactions and indications are encoded using MEDDRA terminology (always?)
# This file attempts to map the reactions in MEDDRA to MESH terms using the CUIs (primarily stored in reaction-meddra.tsv)
# In some cases the term may exist in MEDDRA but not in MESH (hence no cui), so it will store them in unmatched-mesh.tsv
# In those cases where the term was not found in MEDDRA, it will store it in unmatched-reactions.tsv
# Input - mapping file generated by parseFAERS.py (reaction or indication), MESH terminology, MEDDRA terminology in NTriples format, prefix for output files
# Output - three files - reaction_meddra.tsv, unmatched_reactions.tsv, unmatched_mesh.tsv
# python reaction-term.py -i to_use_data_files/faers_reactions.tsv -meddra ../stage_rdf_dump/ontologies/MEDDRA--4 -mesh ../stage_rdf_dump/ontologies/MESH--8 -o to_use_data_files/map -type reaction
# python to_use_scripts/reaction-term.py -i to_use_data_files/faers_indications.tsv -meddra ../stage_rdf_dump/ontologies/MEDDRA--4 -mesh ../stage_rdf_dump/ontologies/MESH--8 -o to_use_data_files/map -type indication

def parse_args():
	p = argparse.ArgumentParser(description='Map the Reactions/Indications in FAERS, generated from parseFAERS.py')
	p.add_argument('-i', '--input', help='Name of the file containing the reactions/indications', default="reaction-faers.tsv")
	p.add_argument('-meddra', '--meddra', help='Location of the N-triple MEDDRA terminology', default="../data/MEDDRA--4")
	p.add_argument('-mesh', '--mesh', help='Location of the N-triple MESH terminology', default="../data/MESH--8")
	p.add_argument('-o', '--prefix', help='prefix of the output file (with directory)', default="results")
	p.add_argument('-type', '--input_type', help='specify if reaction or indication', default="reaction")
	return p.parse_args()

def main():
	args = parse_args()
	meddra_terms = {}
	meddra_ids = {}
	meddrafile = open(args.meddra)
	meddralines = meddrafile.readlines()
	meddrafile.close()

	for k in range(len(meddralines)):
		meddra_parts = meddralines[k].strip().split()
		if "prefLabel" in meddra_parts[1]:
			label = " ".join(meddra_parts[2:])
			label = label.split("@")[0][1:-1]
			print label
			if not label.lower() in meddra_terms:
				meddra_terms[label.lower()] = meddra_parts[0]
		if "cui" in meddra_parts[1]:
			meddra_ids[meddra_parts[0]] = meddra_parts[2].split("^^")[0][1:-1]

	mesh_cuis = {}
	mesh_file = open(args.mesh)
	meshlines = mesh_file.readlines()
	mesh_file.close()

	for k in range(len(meshlines)):
		mesh_parts = meshlines[k].strip().split()
		if "cui" in mesh_parts[1]:
			cui = mesh_parts[2].split("^^")[0][1:-1]
			print cui
			mesh_cuis[cui] = mesh_parts[0][1:-1]

	reaction_meddra = open(args.prefix + "_" + args.input_type + "-meddra.tsv", "w+")
	unmatched_reactions = open(args.prefix + "_unmatched-" + args.input_type + ".tsv", "w+")
	unmatched_mesh = open(args.prefix + "_unmatched-mesh-" + args.input_type + ".tsv", "w+")

	
	def store_matched_reactions(reaction):
		cui = meddra_ids[meddra_terms[reaction.lower()]]
		mesh_id = ""
		if cui in mesh_cuis:
			mesh_id = mesh_cuis[cui]
			reaction_meddra.write(reaction + "\t" + cui + "\t" + mesh_id + "\n")
		else:
			unmatched_mesh.write(reaction + "\t" + cui + "\t" + mesh_id + "\n")

	reactionfile = open(args.input)
	reactionlines = reactionfile.readlines()
	reactionfile.close()

	for r in reactionlines:
		reaction = r.strip()
		if reaction.lower() in meddra_terms:
			store_matched_reactions(reaction)
		elif "^" in reaction:
			reaction = r.strip().replace("^", "'") # strange FAERS struggle
			if reaction.lower() in meddra_terms:
				store_matched_reactions(reaction)
		else:
			unmatched_reactions.write(r)

	unmatched_reactions.close()
	unmatched_mesh.close()
	reaction_meddra.close()

if __name__ == "__main__":
	main()


########################################################################################################################
#                                                                                                                      #
#                      CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data                     #
#                                                                                                                      #
#Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp #
#                                                                                                                      #
#                                              University College London                                               #
#                                                                                                                      #
########################################################################################################################

# This script computes per-residue physicochemical features for the Wuhan-Hu-1 Spike reference sequence:
# hydrophobicity (Kyte-Doolittle scale) net charge at pH 7.4, and molecular weight. Provides helper functions
# for both per-residue values and old -> new substitution deltas (used in sav_features.py to generate in_silico
# features for single AA variant mutations). The main walks the reference Spike sequence and writes a JSON of
# {position -> {charge, molecular_weight, hydrophobicity}}.

from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import json

import numpy as np
from preprocess import filter_seq

kyte_doolittle_scale = {
    'A': 1.8,  # Alanine
    'R': -4.5, # Arginine
    'N': -3.5, # Asparagine
    'D': -3.5, # Aspartic acid
    'C': 2.5,  # Cysteine
    'Q': -3.5, # Glutamine
    'E': -3.5, # Glutamic acid
    'G': -0.4, # Glycine
    'H': -3.2, # Histidine
    'I': 4.5,  # Isoleucine
    'L': 3.8,  # Leucine
    'K': -3.9, # Lysine
    'M': 1.9,  # Methionine
    'F': 2.8,  # Phenylalanine
    'P': -1.6, # Proline
    'S': -0.8, # Serine
    'T': -0.7, # Threonine
    'W': -0.9, # Tryptophan
    'Y': -1.3, # Tyrosine
    'V': 4.2   # Valine
}

pH = 7.4
charge_scores = {}
molecular_weights = {}

for aa in kyte_doolittle_scale.keys():
	analyzed_aa = ProteinAnalysis(aa)
	net_charge = analyzed_aa.charge_at_pH(pH)
	charge_scores[aa] = net_charge
	molecular_weights[aa] = molecular_weight(aa, seq_type='protein')
	
def hydrophobicity_change(old,new):
	return kyte_doolittle_scale[new] - kyte_doolittle_scale[old]

def molecular_weight_change(old,new):
	return molecular_weights[new] - molecular_weights[old]
	
def charge_change(old, new):
	return charge_scores[new] - charge_scores[old]
	
def hydrophobicity_aa(aa):
	return kyte_doolittle_scale[aa]
	
def charge_aa(aa):
	return charge_scores[aa]
	
def molecular_weight_aa(aa):
	return molecular_weights[aa]
	
def reference_values(file_path):
	ref_vals = {}
	with open(file_path, 'r') as f:
		seq = f.readline()
		for i,s in enumerate(seq):
			c = charge_aa(s)
			mw = molecular_weight_aa(s)
			hpb = hydrophobicity_aa(s)
			ref = {"charge":c,"molecular_weight":mw,"hydrophobicity":hpb}
			ref_vals[str(i+1)]= ref
	return ref_vals
	
def mutation_type(mutation):
	if "ins" in mutation or "del" in mutation:
		return "indel"
	else:
		return "sub"
	
def accept_mutation(mutation):
	if mutation_type(mutation) == "sub" and not any(substr in mutation for substr in filter_seq):
		return True
	else:
		return False	
	
def in_silico_matrix(matrix_prefix):
	with open(matrix_prefix + '_mutation_idx.json','r') as f:
		mutation_idx = json.load(f)
	rows = []
	columns = ["hydrophobicity_change","molecular_weight_change","charge_change"]
	for mutation in mutation_idx.keys():
		row = np.full(len(columns),0)
		col_idx = 0
		if accept_mutation(mutation):
			old = mutation[0]
			new = mutation[-1]
			for col in columns:
				func = globals()[col]
				row[col_idx] = func(old,new)
				col_idx += 1
		rows.append(row)	
	
	ins_matrix = np.matrix(rows)
	return columns,ins_matrix

def main():
	ref_vals = reference_values('data/reference/Wuhan-Hu-1-Spike-tidy.aa')
	with open('data/in_silico/Wuhan-Hu-1-Spike.values','w') as f:
		json.dump(ref_vals,f)

if __name__ == '__main__':
    main()

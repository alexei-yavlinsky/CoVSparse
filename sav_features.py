########################################################################################################################
#                                                                                                                      #
#                      CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data                     #
#                                                                                                                      #
#Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp #
#                                                                                                                      #
#                                              University College London                                               #
#                                                                                                                      #
########################################################################################################################

# This script generates per-mutation feature tables for each SARS-CoV-2 protein, used as input to the
# XGBoost model in predict.py. For every accepted single amino-acid variant (excluding insertions,
# deletions, stop codons and disallowed characters), it assembles a row of features combining:
#   - observed mutation count from the GISAID sparse matrix
#   - BLOSUM62 substitution score
#   - in-silico physicochemical changes (charge, hydrophobicity, molecular weight; signed and absolute)
#   - HyPhy selection signatures (class, alpha, beta) [Spike only]
#   - relative solvent accessibility (RSA) [Spike only]
#   - per-site dN/dS [Spike only]
#   - deep mutational scanning phenotypes (ACE2 binding, RBD expression, antibody escape) [Spike only]
# Features that are not available for a given protein are filled with dummy values (typically -100
# or -200). The output is written to data/features/<protein>.csv, one row per mutation.

# The code below was inspired by the following publication:
# Intrahost dynamics, together with genetic and phenotypic effects predict the success of viral mutations
# Cedric CS Tan, Marina Escalera-Zamudio, Alexei Yavlinsky, Lucy van Dorp, Francois Balloux
# bioRxiv 2024.10.18.619070; doi: https://doi.org/10.1101/2024.10.18.619070

from scipy.sparse import csr_matrix, coo_matrix, save_npz, load_npz
import numpy as np
import pandas as pd
from preprocess import proteins, filter_seq
import json
import blosum as bl
from in_silico import charge_change, hydrophobicity_change, molecular_weight_change
from dms import phenotypes
from hyphy import load_standard_definitions
from structure import load_rsa
from dnds import load_dnds
import math
import os

# Count the frequency of each mutation.
# This is the only place we load the raw sparse mutation matrix for single AA variant analysis
# Returns a dictionary mapping mutation names to their total counts in the sparse matrix
# In this and other functions below, matrix_prefix specifies the path to the sparse matrix for a particular protein (e.g. Spike)

def mutation_counts(matrix_prefix):
	matrix = load_npz(matrix_prefix + '.npz') #load the sparse matrix for the given protein
	with open(matrix_prefix + '_mutation_idx.json','r') as f: #get the mappings between sparse matrix column indices and mutation names
		mutation_idx = json.load(f)
	column_names = {value:key for key,value in mutation_idx.items()}
	column_sums = matrix.sum(axis=0) #performing a column-wise sum of the sparse matrix gives us mutation counts
	mutation_counts = {column_names[i]:int(column_sums[0,i]) for i in range(0,column_sums.shape[1])} #create dictionary mapping mutation counts onto mutation names
	return mutation_counts


# Filter function to only accept substitution mutations (exclude insertions, deletions and mutations with disallowed AA characters and stop codons)
def accept_mutation(mutation):
	if mutation_type(mutation) == "sub" and not any(substr in mutation for substr in filter_seq):
		return True
	else:
		return False

# Derive mutation type from AA substitution string
def mutation_type(mutation):
	if "ins" in mutation or "del" in mutation:
		return "indel"
	else:
		return "sub"

# Returns a set of dictionaries, each mapping mutation names to one of hyphy alpha, beta and class values
# This only works for the Spike protein for now		
def hyphy_scores(matrix_prefix,protein):
	default = False
	if protein != 'Spike':
		default = True #if we are not working with the Spike protein, assign default (dummy) values for every mutation
	hyphy_scores = {'hyphy_class':{},'hyphy_alpha':{},'hyphy_beta':{}} #three dictionaries mapping mutation names to hyphy class, alpha and beta values
	definitions = load_standard_definitions(protein) #Spike-only for now, this method will later need to be modified so it accepts multiple proteins
	with open(matrix_prefix + '_mutation_idx.json','r') as f: #load the sparse matrix for the given protein
		mutation_idx = json.load(f)	
	for mutation in mutation_idx.keys():
		if accept_mutation(mutation):			
			class_ = None
			alpha_ = None
			beta_ = None
			if default: #dummy values
				class_ = -200
				alpha_ = -200
				beta_ = -200
			else:
				pos = int(mutation[1:-1]) #extract mutation position
				if pos in definitions: #if this mutation position is covered by hyphy
					def_ = definitions[pos] #load the hyphy definition array for this position
					if def_[0] == 'Diversifying':
						class_ = 1
					elif def_[0] == 'Purifying':
						class_ = -1
					elif def_[0] == 'Neutral':
						class_ = 0
					elif def_[0] == 'Invariant':
						class_ = -100
						print("Invariant dN/dS found")
					alpha_ = def_[2]
					beta_ = def_[3]
				else: #if the mutation position is not covered by hyphy, populate with dummy values
					class_ = -200
					alpha_ = -200
					beta_ = -200
					#print(mutation) #debugging mutation positions not covered by hyphy by printing them to the terminal
			# Populate the entries for this mutation in the hyphy dictionaries
			hyphy_scores['hyphy_class'][mutation] = class_
			hyphy_scores['hyphy_alpha'][mutation] = alpha_
			hyphy_scores['hyphy_beta'][mutation] = beta_
	return hyphy_scores

# Returns a dictionary mapping mutation names to BLOSUM62 scores	
def blosum_scores(matrix_prefix):
	blosum_scores = {}
	blosum62 = bl.BLOSUM(62)
	with open(matrix_prefix + '_mutation_idx.json','r') as f:
		mutation_idx = json.load(f)
	for mutation in mutation_idx.keys():
		if accept_mutation(mutation):
			old = mutation[0]
			new = mutation[-1]
			blosum_scores[mutation] = blosum62[old][new]
	return blosum_scores

# Returns a set of dictionaries, each mapping mutation names onto delta scores for a particular in-silico phenotype
def in_silico_scores(matrix_prefix):
	in_silico_scores = {"charge":{}, "hydrophobicity":{}, "molecular weight":{}, "charge abs":{}, "hydrophobicity abs":{}, "molecular weight abs":{}, }
	with open(matrix_prefix + '_mutation_idx.json','r') as f:
		mutation_idx = json.load(f)
	for mutation in mutation_idx.keys():
		if accept_mutation(mutation):
			old = mutation[0]
			new = mutation[-1]
			in_silico_scores["charge"][mutation] = charge_change(old,new)
			in_silico_scores["hydrophobicity"][mutation] = hydrophobicity_change(old,new)
			in_silico_scores["molecular weight"][mutation] = molecular_weight_change(old,new)
			in_silico_scores["charge abs"][mutation] = abs(charge_change(old,new))
			in_silico_scores["hydrophobicity abs"][mutation] = abs(hydrophobicity_change(old,new))
			in_silico_scores["molecular weight abs"][mutation] = abs(molecular_weight_change(old,new))
	return in_silico_scores

# Returns a dictionary mapping mutation names onto RSA scores, currently for the Spike protein only
def rsa_scores(matrix_prefix,protein):
	default = False
	if protein != 'Spike':
		default = True #if we are not working with the Spike protein, assign default (dummy) values for every mutation
	rsa_scores = {}
	rsa_values = load_rsa() #Spike-only for now, later will need to pass protein to it
	with open(matrix_prefix + '_mutation_idx.json','r') as f:
		mutation_idx = json.load(f)	
	for mutation in mutation_idx.keys():
		if accept_mutation(mutation):
			
			rsa = None
			if default:
				rsa = -100
			else:
				pos = int(mutation[1:-1])
				if pos in rsa_values:
					rsa = rsa_values[pos]
				else:
					rsa = -100
			rsa_scores[mutation] = rsa	
	return rsa_scores

# Returns a dictionary mapping mutation names onto dN/dS values, currently for the Spike protein only	
def dnds_scores(matrix_prefix,protein):
	default = False
	if protein != 'Spike':
		default = True	#if we are not working with the Spike protein, assign default (dummy) values for every mutation
	else:
		syn_values,nonsyn_values,dnds_values = load_dnds(protein)
		inf_val = max([x for x in dnds_values.values() if math.isfinite(x)])*2 #replace infinity values (nonzero divided by zero) with double the maximum finite dN/dS value
	dnds_scores = {}
	nan_val = -100 #replace NaN values (zero divided by zero) with dummy value
	with open(matrix_prefix + '_mutation_idx.json','r') as f:
		mutation_idx = json.load(f)	
	for mutation in mutation_idx.keys():
		if accept_mutation(mutation):
			dnds_value = None
			if default:
				dnds_value = nan_val
			else:
				pos = mutation[1:-1] #extract mutation position
				if pos in dnds_values:
					if syn_values[pos] == 0.0:
						if nonsyn_values[pos] == 0.0:
							dnds_value = nan_val
						else:
							dnds_value = inf_val
					else:
						dnds_value = dnds_values[pos]
			dnds_scores[mutation] = dnds_value
	return dnds_scores

# Returns a set of dictionaries mapping mutation names onto deep mutational scanning (binding, expression) and immune escape phenotypes 	
def phenotype_scores(matrix_prefix,protein):
	phenotype_scores = {}
	default = False
	if protein != 'Spike':
		default = True
	phenotype_names,phenotype_matrix = phenotypes(matrix_prefix)
	num_phenotypes = len(phenotype_names)
	for name in phenotype_names:
		phenotype_scores[name] = {}
	with open(matrix_prefix + '_mutation_idx.json','r') as f:
		mutation_idx = json.load(f)
	for mutation in mutation_idx.keys():
		row = mutation_idx[mutation]
		for i in range(0,num_phenotypes):
			if default:
				value = -100
			else:
				value = phenotype_matrix[row,i]
			phenotype_scores[phenotype_names[i]][mutation] = value
	return phenotype_scores
	

# Generate single AA variant features for each protein and save it into 'data/features/'+protein+'.csv'
def main():
	for protein in proteins:
		print(protein)
		phenotypes_s = phenotype_scores('data/'+protein,protein)
		counts = mutation_counts('data/'+protein)
		blosum = blosum_scores('data/'+protein)
		in_silico = in_silico_scores('data/'+protein)
		hyphy_s = hyphy_scores('data/'+protein,protein)
		rsa_s = rsa_scores('data/'+protein,protein)
		dnds_s = dnds_scores('data/'+protein,protein)
		
		#Create a dictionary of dictionaries, this will be turned into a data frame where each dictionary will be transformed into a column, and its contents into rows
		data = {'count':counts,
				'blosum62':blosum,
				'rsa':rsa_s,
				'dnds':dnds_s,
			   }
		data.update(hyphy_s)
		data.update(phenotypes_s)
		data.update(in_silico)
		df = pd.DataFrame(data)
		filtered_df = df[df.index.map(lambda x: accept_mutation(x))]
		os.makedirs('data/features/', exist_ok=True)
		filtered_df.to_csv('data/features/'+protein+'.csv',index=True)
				
	
if __name__ == '__main__':
    main()

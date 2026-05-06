########################################################################################################################
#                                                                                                                      #
#                      CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data                     #
#                                                                                                                      #
#Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp #
#                                                                                                                      #
#                                              University College London                                               #
#                                                                                                                      #
########################################################################################################################

# Preprocesses deep mutational scanning (DMS) data into per-mutation phenotype features for the SARS-CoV-2 pipeline.
# Reads two source CSVs — Bloom-lab variant binding/expression scores  and antibody-escape fractions — and emits
# per-field JSON files keyed by mutation string. A phenotypes helper then aligns selected fields (escape fraction,
# ACE2 binding, RBD expression) against an external mutation index to build a dense feature matrix with a predefined
# value for missing mutations consumed by the sav_features.py that produces the data/features/{protein}.csv files used
# by predict.py.

# Data sources:
#
# Antibody-escape data (data/dms/escape_fracs.csv):
# Greaney AJ et al. Complete Mapping of Mutations to the SARS-CoV-2 Spike
# Receptor-Binding Domain that Escape Antibody Recognition.
# Cell Host & Microbe. 2021;29(1):44-57.e9.
# doi:10.1016/j.chom.2020.11.007
#
# ACE2-binding and RBD-expression data (data/dms/final_variant_scores.csv):
# Starr TN et al. Deep mutational scans for ACE2 binding, RBD expression, and antibody
# escape in the SARS-CoV-2 Omicron BA.1 and BA.2 receptor-binding domains.
# PLOS Pathogens. 2022;18(11):e1010951.
# doi:10.1371/journal.ppat.1010951

import numpy as np
import pandas as pd
import json

variants = ["Alpha","Beta","Delta","Eta","Omicron_BA1","Omicron_BA2","Wuhan-Hu-1_v1","Wuhan-Hu-1_v2"]
fields = ["bind","delta_bind","n_bc_bind","n_libs_bind","bind_rep1","bind_rep2","bind_rep3","expr","delta_expr","n_bc_expr","n_libs_expr","expr_rep1","expr_rep2"]
field_subset = ["bind","expr"] #delta_bind,delta_expr

esc_fields = ["mut_escape_frac_epistasis_model","mut_escape_frac_single_mut","site_total_escape_frac_epistasis_model","site_total_escape_frac_single_mut","site_avg_escape_frac_epistasis_model","site_avg_escape_frac_single_mut"]
esc_field_subset = ["mut_escape_frac_single_mut"] #mut_escape_frac_epistasis_model

monoclonals = ["COV2-2050_400","COV2-2082_400","COV2-2094_400","COV2-2096_400","COV2-2165_400","COV2-2479_400","COV2-2499_400","COV2-2677_400","COV2-2832_400","CR3022_400"]

def append_mutation(row):
	return row["wildtype"] + str(row["protein_site"]) + row["mutation"]
	

def binding():
	bind = pd.read_csv("data/dms/final_variant_scores.csv",low_memory=False)
	for variant in variants:
		_bind = bind[bind['target']==variant]
		for field in fields:
			js = _bind[["mutation",field]].to_json()
			with open("data/dms/bind_"+variant+"_"+field+".json", "w") as write_file:
				write_file.write(js)	

def escape():
	agg_dict = {}
	for field in esc_fields:
		agg_dict[field] = 'sum'
	agg_dict['selection'] = 'count'
	
	escape_vals = pd.read_csv("data/dms/escape_fracs.csv",low_memory=False)
	escape_vals = escape_vals[escape_vals["library"] == "average"]
	print(escape_vals.shape[0])
	escape_vals['mut_string'] = escape_vals.apply(lambda row: append_mutation(row), axis = 1)
	agg = escape_vals.groupby("mut_string").agg(agg_dict).reset_index() 
	for field in esc_fields:
		js = agg[["mut_string",field]].to_json()
		with open("data/dms/esc_"+field+".json", "w") as write_file:
			write_file.write(js)
			
def reorder_fill_row(mutation_idx, mutation_values,field,na_val=-100):
	#print("-----")
	#row = np.zeros(len(mutation_idx))
	row = np.full(len(mutation_idx),na_val)
	for key in mutation_values['mutation'].keys():
		mutation = mutation_values['mutation'][key]
		value = mutation_values[field][key]
		col_index = mutation_idx.get(mutation)
		if col_index is None:
			#print(mutation)
			continue
		if value is not None:
			row[col_index] = value
	return row
	
def phenotypes(matrix_prefix,na_val_esc=-1,na_val_bind=-100):
	variant = "Wuhan-Hu-1_v1"
	with open(matrix_prefix + '_mutation_idx.json','r') as f:
		mutation_idx = json.load(f)
	rows = []
	columns = []
	for field in esc_field_subset:
		with open("data/dms/esc_"+field+".json",'r') as f:
			mutation_values = json.load(f)
			mutation_values['mutation'] = mutation_values.pop('mut_string')
			rows.append(reorder_fill_row(mutation_idx,mutation_values,field,na_val=na_val_esc))
			columns.append(field)
			
	for field in field_subset:
		with open("data/dms/bind_"+variant+"_"+field+".json",'r') as f:
			mutation_values = json.load(f)
			rows.append(reorder_fill_row(mutation_idx,mutation_values,field,na_val=na_val_bind))
			columns.append(field)
	phenotypes = np.matrix(rows).transpose()		
	#phenotypes_csr = csr_matrix(phenotypes)		
	return columns,phenotypes
				
def main():
	binding()
	escape()

if __name__ == '__main__':
    main()
    
   #bloom et al fitness effects data

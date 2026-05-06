########################################################################################################################
#                                                                                                                      #
#                      CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data                     #
#                                                                                                                      #
#Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp #
#                                                                                                                      #
#                                              University College London                                               #
#                                                                                                                      #
########################################################################################################################

# This script constructs a sparse matrix for each SARS-CoV-2 protein

from preprocess import (
	proteins,
	get_columns,
	get_substitution_col,
	get_accession_col,
	get_clade_col,
	get_date_col,
	get_variant_col,
	get_location_col,
)
from scipy.sparse import csr_matrix, coo_matrix, save_npz, load_npz
import numpy as np
import gzip
from functools import cmp_to_key
import re
import json

# Helper method for sorting mutations by position
# Strips non-numeric characters with a regex to extract the position number to be used with cmp_to_key for custom sorting 
def compare_mutations(a,b):
	a_num = int(re.sub("[^0-9]","",a))
	b_num = int(re.sub("[^0-9]","",b))
	if a_num > b_num:
		return 1
	elif a_num < b_num:
		return -1
	else:
		return 0
		
		
#Build a sparse matrix from a protein-specific metadata file at input_path
#Saves the sparse matrix and other associated files to output_prefix
def build_matrix(input_path,output_prefix):
	mutations = set([])
	mutation_idx = {}
	substitution_col = get_substitution_col()
	accession_col = get_accession_col()
	clade_col = get_clade_col()
	date_col = get_date_col()
	variant_col = get_variant_col()
	location_col = get_location_col()
	columns = get_columns()
	substitution_col_idx = columns.index(substitution_col)
	accession_col_idx = columns.index(accession_col)
	clade_col_idx = columns.index(clade_col)
	date_col_idx = columns.index(date_col)
	location_col_idx = columns.index(location_col)
	variant_col_idx = columns.index(variant_col)
	
	f = gzip.open(input_path,'rt')
	line_counter = 0
	
	n_nonzero = 0 # counter for the number of nonzero elements
	entries = [] #array for sequence metadata values 
	
	# First, scan the metadata file to create the list of all unique mutations,
	# count the number of nonzero elements and store other metadata fields
	for line in f:
		if (line_counter == 0):
			pass #skip the header
		else:
			values = line.strip().split('\t')
			mutation_list = values[substitution_col_idx].split(',')
			mutation_set = set(mutation_list)
			if '' in mutation_set:
				mutation_set.remove('')
			mutations.update(mutation_set)
			n_nonzero += len(mutation_set) #increment the number of nonzero elements by the number of mutations encountered for this sequence
			entries.append([values[accession_col_idx],values[clade_col_idx],values[variant_col_idx],values[date_col_idx],values[location_col_idx]]) #store other metadata for this sequence
		line_counter = line_counter + 1
	f.close()
	unique_mutations = list(mutations)
	mutations_cmp_key = cmp_to_key(compare_mutations)
	unique_mutations.sort(key=mutations_cmp_key) #sort mutations by position
	
	# The scipy sparse matrix object used in this implementation does not have column names, only indices
	# Below, we are creating a table linking column names (mutations) with column indices
	idx = 0
	for um in unique_mutations:
		mutation_idx[um] = idx
		idx = idx + 1
	  
	n_items = len(entries)
	n_mutations = len(unique_mutations)
	matrix_stats = {"Non-zero entries": n_nonzero, "Number of entries": n_items, "Unique mutations":n_mutations}
	print(matrix_stats)
	
	# Create the data structures and allocate storage for the sparse matrix, using the number of nonzero elements
	data = np.empty(n_nonzero,dtype=np.intc) #nonzero mutation counts
	rows = np.empty(n_nonzero,dtype=np.intc) #row indices for the mutation coutns
	cols = np.empty(n_nonzero,dtype=np.intc) #column indices for the mutation counts
	
	# Read the metadata file again, this time populating the sparse matrix as we go along
	f = gzip.open(input_path,'rt')
	ind = 0 #index counter for where we are in the sparse matrix, this is used for populating the sparse matrix data structures slice by slice (where each slice are the nonzero values of a row === nonzero mutation counts in the sequence)
	line_counter = 0 #row counter
	for line in f:
		if (line_counter % 100000 == 0):
			print(line_counter)
		if (line_counter == 0):
			pass #skip the header (as in the first pass)
		else:
			values = line.strip().split('\t')
			mutation_list = values[substitution_col_idx].split(',')
			mutation_set = set(mutation_list)
			if '' in mutation_set:
				mutation_set.remove('')
			mutation_indices = []
			# Look up column indices for the mutations encountered in this sequence
			for m in mutation_set:
				mutation_indices.append(mutation_idx[m])
			mutation_counts = [1 for x in mutation_indices] #a mutation can only occur once in a sequence, so all counts are set to 1
			mutation_entries = [line_counter-1 for x in mutation_indices] #for all mutations in this sequence, set its row index as the current row number
			n_vals = len(mutation_indices) #the number of non-zero entries for this sequence/row
			ind_end = ind + n_vals #index counter of where we will be after we process this sequence/row
			data[ind:ind_end] = mutation_counts #populate the portion of the nonzero mutation count data structure corresponding to this row  
			cols[ind:ind_end] = mutation_indices #write the column indices for the nonzero mutation count data stored above
			rows[ind:ind_end] = mutation_entries #write the row indices for the nonzero mutatation count data
			ind = ind_end #shift the sparse matrix index counter to the end of the data for the current row
		line_counter = line_counter + 1
	f.close()
	# Create scipy's coo_matrix using the above sparse matrix data structures and then convert it to csr_matrix and save it
	cm = coo_matrix((data, (rows,cols)), shape=(n_items,n_mutations),dtype=np.intc)
	csr = cm.tocsr()
	save_npz(output_prefix + '.npz',csr)
	
	# Save the corresponding metadata, mappings between column indices and mutation names, and summary statistics for this matrix
	with gzip.open(output_prefix + '_entries.json.gz','wt') as f:
		json.dump(entries,f)
	with open(output_prefix + '_mutation_idx.json','w') as f:
		json.dump(mutation_idx,f)
	with open(output_prefix + '_stats.json','w') as f:
		json.dump(matrix_stats,f)

# Build a sparse matrix for each protein
def build_matrices():
	for protein in proteins:
		print(protein)
		build_matrix('data/'+protein+'.tsv.gz','data/'+protein)
	
def main():
	build_matrices()
	
if __name__ == '__main__':
    main()

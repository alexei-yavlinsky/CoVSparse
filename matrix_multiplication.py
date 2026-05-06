########################################################################################################################
#                                                                                                                      #
#                      CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data                     #
#                                                                                                                      #
#Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp #
#                                                                                                                      #
#                                              University College London                                               #
#                                                                                                                      #
########################################################################################################################

# Demonstrates how sparse matrix multiplication can aggregate per-mutation scores up to per-sequence summaries
# across the whole GISAID dataset in one shot. For a given sequence x mutation sparse matrix, it multiplies by
# a mutation x phenotype matrix (DMS binding/expression/escape) and a mutation x in-silico matrix
# (charge/hydrophobicity/molecular-weight deltas), producing a sequence x feature matrix where each cell is the
# sum of the relevant per-mutation scores across all mutations carried by that sequence.

from scipy.sparse import csr_matrix, coo_matrix, save_npz, load_npz, hstack
import numpy as np
import json
from in_silico import in_silico_matrix
from dms import phenotypes


def additive_features(matrix_prefix):
	matrix = load_npz(matrix_prefix + '.npz') #load the sparse matrix for the given protein
	
	phenotype_names,phenotype_matrix = phenotypes(matrix_prefix,na_val_esc=0,na_val_bind=0) #get the matrix of phenotype values, where each row is a mutation and each column is a DMS phenotype
	phenotypes_csr_matrix = csr_matrix(phenotype_matrix) #turn it into a sparse matrix
	phenotypes_product = matrix * phenotypes_csr_matrix #perform dot product between the GISAID sequence sparse matrix and the phenotype value matrix
	print(phenotypes_product.shape)
	
	ins_names,ins_matrix = in_silico_matrix(matrix_prefix) #get the matrix of in-silico values, where each row is a mutation and each column is the difference in in-silico scores between the old and new AA
	ins_csr_matrix = csr_matrix(ins_matrix)  #turn it into a sparse matrix
	ins_product = matrix * ins_csr_matrix #perform dot product between the GISAID sequence sparse matrix and the in-silico value matrix
	print(ins_product.shape)
	
	full_matrix = hstack([phenotypes_product,ins_product],format="csr")
	print(full_matrix.shape)
	
	return full_matrix

def main():
	full_matrix = additive_features('data/Spike')
	
if __name__ == '__main__':
    main()

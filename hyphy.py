########################################################################################################################
#                                                                                                                      #
#                      CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data                     #
#                                                                                                                      #
#Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp #
#                                                                                                                      #
#                                              University College London                                               #
#                                                                                                                      #
########################################################################################################################

# This script loads data derived from site-level FEL selection-analysis output
# for the SARS-CoV-2 Spike protein using HyPhy (Kosakovsky Pond et al.), distributed at
# https://github.com/spond/SARS-CoV-2-variation

import sys
import tarfile
import lzma
import gzip
import csv

def load_definitions(filename):
	definitions = {}
	with open(filename, 'r') as f:
		reader = csv.DictReader(f, delimiter=',')
		for row in reader:
			site_ = row['Codon'].split('/')[1]
			class_ = row['Class']
			pval_ = row['p-value']
			alpha_ = row['α']
			beta_ = row['β']
			definitions[int(site_)] = (class_,float(pval_),alpha_,beta_)
			#print(site_,class_,pval_)
	return definitions		

def load_standard_definitions(protein):
	return load_definitions('data/hyphy/Spike_dNdS_07_2023.csv')
							
def main():
	definitions = load_definitions('data/hyphy/Spike_dNdS_07_2023.csv')
	print(definitions)

if __name__ == '__main__':
    main()

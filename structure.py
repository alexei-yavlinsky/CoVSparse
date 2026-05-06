########################################################################################################################
#                                                                                                                      #
#                      CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data                     #
#                                                                                                                      #
#Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp #
#                                                                                                                      #
#                                              University College London                                               #
#                                                                                                                      #
########################################################################################################################

# This script loads the RSA data from data/structure/6vxx.filt.SASA.csv: per-residue solvent-accessible surface
# area for the SARS-CoV-2 Spike trimer in the closed conformation, computed from PDB entry 6VXX
# (Walls et al. 2020, Cell 181:281-292; doi:10.1016/j.cell.2020.02.058).

import sys
import tarfile
import lzma
import gzip
import csv

def load_rsa(): #s
	definitions = {}
	rsa_values = {}
	with open('data/structure/6vxx.filt.SASA.csv', 'r') as f:
		reader = csv.DictReader(f, delimiter=',')
		for row in reader:
			site_ = int(row['codon_number'])
			chain_ = row['chain']
			sasa_ = row['SASA']
			maxAsa_ = row['maxASA']
			if site_ not in definitions:
				definitions[site_] = {}
			definitions[site_][chain_] = float(sasa_)/float(maxAsa_)
	
	for site_ in definitions.keys():
		chains = definitions[site_]
		rsa = 0
		for chain_ in chains.keys():
			rsa += chains[chain_]
		rsa = rsa/len(chains)
		rsa_values[site_] = rsa	
	return rsa_values


def main():
	rsa_values = load_rsa()
	print(rsa_values)

if __name__ == '__main__':
    main()

########################################################################################################################
#                                                                                                                      #
#                      CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data                     #
#                                                                                                                      #
#Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp #
#                                                                                                                      #
#                                              University College London                                               #
#                                                                                                                      #
########################################################################################################################

# This script loads the dN/dS data (computed externally) for each site in the corresponding SARS-CoV-2 protein

import json		

def load_dnds(protein):
	with open(f'data/dnds/{protein}_syn.json','r') as syn, open(f'data/dnds/{protein}_nonsyn.json','r') as nonsyn, open(f'data/dnds/{protein}_dnds.json','r') as ratio:
		syn_dict = json.load(syn)
		nonsyn_dict = json.load(nonsyn)
		dnds_dict = json.load(ratio)
		return syn_dict, nonsyn_dict, dnds_dict
							
def main():
	pass

if __name__ == '__main__':
    main()

########################################################################################################################
#                                                                                                                      #
#                      CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data                     #
#                                                                                                                      #
#Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp #
#                                                                                                                      #
#                                              University College London                                               #
#                                                                                                                      #
########################################################################################################################

# This script filters GISAID SARS-CoV-2 metadata and prepares it for sparse matrix construction

import sys
import tarfile
import lzma
import gzip

#Columns to keep from the metadata file
cols_to_keep = ["Accession ID","Collection date","Location","Host","Clade","Pango lineage","Variant","AA Substitutions","Is complete?","Is high coverage?","Is low coverage?"]

#Which sequences to consider
filter_values = {"Host":"Human","Is low coverage?":"","Is complete?":"True"}

#Column names in the metadata file
substitution_col = "AA Substitutions"
accession_col = "Accession ID"
clade_col = "Clade"
date_col = "Collection date"
lineage_col = "Pango lineage"
variant_col = "Variant"
location_col = "Location"

#Proteins to be considered
proteins = ["E","M","N","NS3","NS6","NS7a","NS7b","NS8","NSP1","NSP10","NSP11","NSP12","NSP13","NSP14","NSP15","NSP16","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","Spike"]
protein_set = set(proteins)

#Ignore mutations that have the following AA values
filter_seq = ['B','Z','J','X','O','stop'] 

#Protein regions (Spike only for now)
protein_regions =	{
						"Spike" : 	{
										"SARS-CoV-like_Spike_S1_NTD": [13,304],
										"SARS-CoV-2_Spike_S1_RBD": [319,541],
										"SARS-CoV-like_Spike_SD1-2_S1-S2_S2": [543,1208]
									}
					}

#Helper methods to get names of appropriate columns in the metadata file

def get_columns():
	return cols_to_keep
	
def get_substitution_col():
	return substitution_col
	
def get_accession_col():
	return accession_col
	
def get_clade_col():
	return clade_col
	
def get_date_col():
	return date_col
	
def get_lineage_col():
	return lineage_col
	
def get_variant_col():
	return variant_col
	
def get_location_col():
	return location_col

# Strip out the protein prefix (e.g. Spike from Spike_L455S) and insert
# the substitution to the corresponding protein array

def extract_substitutions(value):
	proc_value = value.lstrip('(').rstrip(')')
	proc_tokens = proc_value.split(',')
	filtered_tokens = {protein: [] for protein in proteins}
	for token in proc_tokens:
		if token == '':
			continue
		prot_sub = token.split('_')
		protein = prot_sub[0]
		sub = prot_sub[1]
		if protein in protein_set:
			filtered_tokens[protein].append(sub)
	return filtered_tokens
	
# Method for processing each sequence line in the metadata file

def filter_line(line,col_indices):
	values = line.strip().split('\t')
	if values[col_indices["Clade"]] == '':
		return None #if we don't have clade information, ignore the sequence
	extracted_values = {protein: [] for protein in proteins}
	for col in cols_to_keep:
		value = values[col_indices[col]]
		#if we are filtering on column col and its value doesn't match our criteria, ignore the sequence
		if col in filter_values and value != filter_values[col]:
			return None
		#extract substitution values and place them in appropriate protein dictionary entries
		if col == substitution_col:
			sub_values = extract_substitutions(value)
			for protein in proteins:
				extracted_values[protein].append(','.join(sub_values[protein]))
		else:
			for protein in proteins: 
				extracted_values[protein].append(value)
	output_lines = {protein:'\t'.join(extracted_values[protein]) for protein in proteins}
	
	#for each input sequence, return a dictionary where,
	#for each protein, we have substitutions specific to that protein
	#plus a copy of all the other metadata
	#
	#thus, for each input line in the GISAID metadata file, we have
	#one output line, per protein, that will be stored in the corresponding
	#protein-specific output file
	
	return output_lines 
	
# Read the column names from the GISAID metadata file to get their positional indices	
def read_header(header):
	column_names = header.strip().split('\t')
	col_indices = {}
	for i in range(0,len(column_names)):
		if column_names[i] in cols_to_keep:
			col_indices[column_names[i]] = i
	return col_indices	

# Main method to read the GISAID metadata file
# metadata_only=True --> ignore the AA substitution data and only print the sequence ID to the terminal
def process_gisaid(gisaid_path,gisaid_output_path,metadata_only=False): 
	col_indices = None
	line_counter = 0
	#set up the output metadata tsv files, one file per protein
	if not metadata_only:
		w = {protein: gzip.open(gisaid_output_path+'/'+protein+'.tsv.gz','wt') for protein in proteins}
	with lzma.open(gisaid_path, mode="rb") as lzma_file:
		with tarfile.open(fileobj=lzma_file, mode="r|") as tar:
			for member in tar:
				print(f"Found file: {member.name}",file=sys.stderr)
				if member.name == "metadata.tsv":
					f = tar.extractfile(member)
					if f is not None:
						for binary_line in f: #process each line in the input metadata file
							line = binary_line.decode('utf-8')
							if (line_counter % 1000 == 0):
								print(line_counter,file=sys.stderr)
							if (line_counter == 0):
								col_indices = read_header(line)
								if not metadata_only:
									for protein in proteins: #write the column names in each protein tsv file
										w[protein].write('\t'.join(cols_to_keep)+'\n')
							else:
								output_lines = filter_line(line,col_indices)
								if output_lines is not None:
									if not metadata_only:
										for protein in proteins: #write the protein-specific metadata for the input sequence into the corresponding protein tsv file
											w[protein].write(output_lines[protein]+'\n')
									else:
										first_protein = output_lines[list(output_lines.keys())[0]]
										sequence_id = first_protein.split("\t")[0]
										print(sequence_id) 
							line_counter = line_counter + 1
	if not metadata_only:
		for protein in proteins:
			w[protein].close()
	
def main():
	process_gisaid('data/metadata_tsv.tar.xz','data/') #read compressed data directly


if __name__ == '__main__':
    main()

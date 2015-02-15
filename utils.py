from pprint import pprint

fasta_filenames = ["data/HIV-1_gag.fasta",
				   "data/HIV-1_nef.fasta"]

def read_fasta_file(fasta_file, aligned = True):
	sequences = []
	curr_seq = ""
	for line in open(fasta_file, 'r').readlines():
		if '>' in line:
			if len(curr_seq) != 0:
				sequences.append(curr_seq)
				curr_seq = ''
		else:
			curr_seq += line.replace('\n', '').strip()

	# Remove "-" placeholders that code for empty spaces in the aligned
	# FASTA file
	if not aligned:
		for i in xrange(len(sequences)):
			sequences[i] = sequences[i].replace('-', '')

	return sequences

def get_all_sequences(aligned = True):
	all_seqs = []
	for file in fasta_filenames:
		all_seqs.extend(read_fasta_file(file, aligned))
	return all_seqs		

from pprint import pprint

fasta_filenames = ["data/HIV-1_gag.fasta",
				   "data/HIV-1_nef.fasta"]

epitope_coverage = {} # epitope : coverage
epitope_length = 9

def coverage(mosaic_seq, population_seqs):
	""" Returns metric for coverage of a mosaic sequence based on the
	    population.

		For right now, this algorithm iterates through sliding windows of 9mers
		for each protein in the datain the mosaic protein and """

	coverage = 0
	already_seen = set() # Avoids double counting coverage (within a mosaic protein)
	for start_i in xrange(len(mosaic_seq) - epitope_length + 1):
		curr_epitope = mosaic_seq[start_i:start_i + epitope_length]
		if curr_epitope in already_seen:
			continue
		else:
			already_seen.add(curr_epitope)
			
		if curr_epitope in epitope_coverage:
			# Coverage is cached
			coverage += epitope_coverage[curr_epitope]
		else:
			# Need to calculate coverage by iterating through all sequences.
			# In each iteration, calculate number of distinct epitopes.
			# If epitope in natural sequence, coverage = 1.0 / number of epitopes in sequence
			epitope_coverage[curr_epitope] = 0.0
			for seq in population_seqs:
				curr_natural_seq_epitopes = set()
				epitope_match = 0.0
				curr_coverage = 0.0
				for seq_start_i in xrange(len(seq) - epitope_length + 1):
					if seq[seq_start_i:seq_start_i + epitope_length] == curr_epitope:
						epitope_match = 1.0
					curr_natural_seq_epitopes.add(seq[seq_start_i:seq_start_i + epitope_length])
				epitope_coverage[curr_epitope] += epitope_match / len(curr_natural_seq_epitopes)
			epitope_coverage[curr_epitope] /= len(population_seqs)
			coverage += epitope_coverage[curr_epitope]
	return coverage

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

	if not aligned:
		for i in xrange(len(sequences)):
			sequences[i] = sequences[i].replace('-', '')

	return sequences

def get_all_sequences(aligned = True):
	""" This function returns a list of all fasta sequences for HIV from
        the global list, fasta_filenames.

		Remove '-' placeholders that code for empty spaces in the FASTA files
    """
	
	all_seqs = []
	for file in fasta_filenames:
		all_seqs.extend(read_fasta_file(file, aligned))
	return all_seqs		

if __name__ == "__main__":
	# Coverage test, should output 0.67
	print coverage('ABCDEFGHIJKL', ['ABCDEFGHIJKLMN'])

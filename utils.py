from pprint import pprint

fasta_filenames = ["data/SHORT_HIV-1_gag.fasta",
                   "data/SHORT_HIV-1_nef.fasta"]

soft_epitope_coverage = {} # epitope : coverage
hard_epitope_coverage = {} # epitope : coverage
epitope_length = 9

def coverage(mosaic_seq, population_seqs, t="soft"):
    """ Returns metric for coverage of a mosaic sequence based on the
        population.

        t (type) parameter: either 'soft' or 'hard' matching for epitopes.  In other words,
        soft matching takes into account how many amino acids match in the same positions
        between two 9mers, while hard matching requires the entire sequence to be equivalent
        to contribute to the coverage score. 'hard' corresponds to Fisher's paper.

        NOTE: the coverage metric doesn't really take into account insertions/deletions; there
        could be an 8mer that's really close to an epitope but doesn't get primed for additio
        b/c our metric would not pick that up as well. """

    coverage = 0
    already_seen = set() # Avoids double counting coverage (within a mosaic protein)
    for start_i in xrange(len(mosaic_seq) - epitope_length + 1):
        curr_epitope = mosaic_seq[start_i:start_i + epitope_length]
        if curr_epitope in already_seen:
            continue
        else:
            already_seen.add(curr_epitope)

        if t == 'soft' and curr_epitope in soft_epitope_coverage:
            # Coverage is cached
            coverage += soft_epitope_coverage[curr_epitope]
        elif t == 'hard' and curr_epitope in hard_epitope_coverage:
            coverage += hard_epitope_coverage[curr_epitope]
        else:
            # Need to calculate coverage by iterating through all sequences.
            # In each iteration, calculate number of distinct epitopes.
            # If epitope in natural sequence, coverage = 1.0 / number of epitopes in sequence
            if t == 'soft':
                soft_epitope_coverage[curr_epitope] = 0.0
            else:
                hard_epitope_coverage[curr_epitope] = 0.0
                
            for seq in population_seqs:
                curr_natural_seq_epitopes = set()
                epitope_match = 0.0
                curr_coverage = 0.0
                for seq_start_i in xrange(len(seq) - epitope_length + 1):
                    if t == 'hard':
                        if seq[seq_start_i:seq_start_i + epitope_length] == curr_epitope:
                            epitope_match = 1.0
                    else: # type == 'soft'
                        aa_position_match_count = 0.0
                        for aa_position in xrange(epitope_length):
                            if curr_epitope[aa_position] == seq[seq_start_i + aa_position]:
                                aa_position_match_count += 1.0
                        if aa_position_match_count / epitope_length > epitope_match:
                            epitope_match = aa_position_match_count / epitope_length
                        
                    curr_natural_seq_epitopes.add(seq[seq_start_i:seq_start_i + epitope_length])
                if t == 'soft':
                    soft_epitope_coverage[curr_epitope] += epitope_match / len(curr_natural_seq_epitopes)
                else:
                    hard_epitope_coverage[curr_epitope] += epitope_match / len(curr_natural_seq_epitopes)

            if t == 'soft':
                soft_epitope_coverage[curr_epitope] /= len(population_seqs)
                coverage += soft_epitope_coverage[curr_epitope]
            else:
                hard_epitope_coverage[curr_epitope] /= len(population_seqs)
                coverage += hard_epitope_coverage[curr_epitope] 
            
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

        Remove '-' placeholders that code for empty spaces in the FASTA files """
    
    all_seqs = []
    for file in fasta_filenames:
        all_seqs.extend(read_fasta_file(file, aligned))
    return all_seqs     

if __name__ == "__main__":
    # Coverage test, should output 0.67 (hard)
    print coverage('ABCDEFGHIJKL', ['ABCDEFGHIJKLMN', 'ABCDEFGHIJKLMN'], 'hard')
    # Coverage test, should be 0.1111 (soft)
    print coverage('ABCDEFGHI', ['ATRYSEFSG'], 'soft')

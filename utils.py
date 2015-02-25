from pprint import pprint
import random

fasta_filenames = ["data/HIV-1_gag.fasta",
                   "data/HIV-1_nef.fasta"]

soft_epitope_coverage = {} # epitope : coverage
hard_epitope_coverage = {} # epitope : coverage
epitope_length = 9

possible_mutations = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

aa_ngrams = {} # prior & after : middle : count
aa_top_ngrams = {} # prior & after : middle


def write_coverage_to_file(filename, cov_type = 'soft'):
    with open(filename, 'w') as f:
        dict_to_write = None
        if cov_type == 'soft':
            dict_to_write = soft_epitope_coverage
        else:
            dict_to_write = hard_epitope_coverage
        for epitope in dict_to_write:
            f.write(epitope + ',' + str(dict_to_write[epitope]) + '\n')
    
def read_coverage_from_file(filename, cov_type = 'soft'):
    lines = open(filename, 'r').readlines()
    for l in lines:
        epitope, coverage = l.split(',')
        if cov_type == 'soft':
            soft_epitope_coverage[epitope] = float(coverage)
        else:
            hard_epitope_coverage[epitope] = float(coverage)
        
def choose_mutation(mosaic_seq, init_coverage, population):
    # Handle beginning and end
    # here. soon.

    # Rest of the amino acids
    top_choices = [] # (position, letter, coverage)
    for i in xrange(1, len(mosaic_seq) - 1):
        # Look up neighbors.  If not there, use default probability of 0.0001
        neighbors = (mosaic_seq[i-1], mosaic_seq[i+1])
        max_mutations = 5
        mutation_choices = None
        if neighbors in aa_top_ngrams:
			num_mutations = max(max_mutations, len(aa_top_ngrams[neighbors]))
			mutation_choices = [aa_top_ngrams[neighbors][m][0] for m in xrange(num_mutations)]
        else:
            mutation_choices = [possible_mutations[int(random.random() * 20)] for m in xrange(max_mutations)]

        print mutation_choices
        for mutation_choice in mutation_choices:
            if mutation_choice != mosaic_seq[i]:
                mutated_sequence = mosaic_seq[:i] + mutation_choice + mosaic_seq[i+1:]
                curr_coverage = coverage(mutated_sequence, population)
                if curr_coverage > init_coverage:
                    top_choices.append((i, mutation_choice, curr_coverage))

    if len(top_choices) == 0:
        return (-1, None, init_coverage)
    else:
        num_top_choices_considered = 3
        top_choices = sorted(top_choices, key=lambda x: x[2], reverse=True)
        num_considered = min(len(top_choices), num_top_choices_considered)
        return top_choices[int(random.random() * num_considered)]

def calc_aa_ngrams(pop_seqs):
    for seq in pop_seqs:
        for i in xrange(len(seq) - 2):
            curr_chunk = seq[i:i+3]
            neighbors = (curr_chunk[0], curr_chunk[2])
            if neighbors not in aa_ngrams:
                aa_ngrams[neighbors] = {}
            if curr_chunk[1] not in aa_ngrams[neighbors]:
                aa_ngrams[neighbors][curr_chunk[1]] = 0.0
            aa_ngrams[neighbors][curr_chunk[1]] += 1.0
            
        beg_neighbor = ("_", seq[1])
        if beg_neighbor not in aa_ngrams:
            aa_ngrams[beg_neighbor] = {}
        if seq[0] not in aa_ngrams[beg_neighbor]:
            aa_ngrams[beg_neighbor][seq[0]] = 0.0
        aa_ngrams[beg_neighbor][seq[0]] += 1.0

        end_neighbor = (seq[-2], "_")
        if end_neighbor not in aa_ngrams:
            aa_ngrams[end_neighbor] = {}
        if seq[-1] not in aa_ngrams[end_neighbor]:
            aa_ngrams[end_neighbor][seq[-1]] = 0.0
        aa_ngrams[end_neighbor][seq[-1]] += 1.0

    for neighbor in aa_ngrams:
        curr_sum = 0.0
        for middle in aa_ngrams[neighbor]:
            curr_sum += aa_ngrams[neighbor][middle]

        for middle in aa_ngrams[neighbor]:
            aa_ngrams[neighbor][middle] /= curr_sum

    for neighbor in aa_ngrams:
        middle_probabilities = [(middle, aa_ngrams[neighbor][middle]) for middle in aa_ngrams[neighbor]]
        aa_top_ngrams[neighbor] = sorted(middle_probabilities, key=lambda x: x[1], reverse=True)[:5]

def random_mutation(sequence):
    position =  int(random.random() * sequence)
    amino_acid = possible_mutations[int(random.random() * 20)]

    return (position, amino_acid)

""" DEPRECATED FUNCTION, stored for backup and/or comparison purposes """
def get_location_probabilities(mosaic_seq):
    """ Scores amino acid by the mean of the coverage of nearby epitopes

    At first, this will choose to point mutate amino acids with no epitopes around it.
    Once most amino acids have epitopes nearby, it will start to mutate problem amino acids
    in viable epitopes.

    Problem: if we can't viably mutate amino acids in one section, then we'll be stuck.
    Solution: either make the choice of mutation location probabilistic OR
    keep counts and factor that count score in as well """
    
    amino_acid_importance_scores = [0.0] * len(mosaic_seq)
    for i in xrange(len(mosaic_seq)):
        sum_sliding_window_epitope_coverage_scores = 0.0
        first_epitope_start_i = max(i - epitope_length + 1, 0)
        last_epitope_start_i = min(i, len(mosaic_seq) - epitope_length)
        for epitope_start_i in xrange(first_epitope_start_i, last_epitope_start_i + 1):
            epitope_end_i = epitope_start_i + epitope_length
            epitope_seq = mosaic_seq[epitope_start_i:epitope_end_i]
            sum_sliding_window_epitope_coverage_scores += soft_epitope_coverage[epitope_seq]
        num_epitopes = last_epitope_start_i - first_epitope_start_i + 1
        amino_acid_importance_scores[i] = sum_sliding_window_epitope_coverage_scores / num_epitopes

    # Choose random point mutation
    amino_acid = possible_mutations[int(random.random() * 20)]  
    location = amino_acid_importance_scores.index(min(amino_acid_importance_scores))

    return location, amino_acid

    

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
                            break
                    else: # type == 'soft'
                        aa_position_match_count = 0.0
                        for aa_position in xrange(epitope_length):
                            if curr_epitope[aa_position] == seq[seq_start_i + aa_position]:
                                aa_position_match_count += 1.0
                        perfect_match_weight = 2.0
                        if aa_position_match_count == epitope_length:
                            epitope_match += perfect_match_weight # To weight perfect epitope matches highery
                        else:
                            epitope_match += aa_position_match_count / epitope_length
                        
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
    #print coverage('ABCDEFGHIJKL', ['ABCDEFGHIJKLMN', 'ABCDEFGHIJKLMN'], 'hard')
    #print coverage('ABCDEFGHIJKLMNOP', ['ATRYSEFSG'], 'soft')

    #print choose_mutation('ABCDEFGHIJKLMNOP')
    #calc_aa_ngrams(['ABC', 'ADC'])
    # pop = get_all_sequences(aligned=False)
    # read_coverage_from_file("test_coverage_file.txt")
    # init_coverage = coverage('TYKGAFDLSHFLKEKGGLDGLVWSPKRQEILDLWVYHTQGYFPDWQNYTP', pop)
    # print "Calculated coverage once, time starting now"
    # calc_aa_ngrams(pop)
    # print choose_mutation('TYKGAFDLSHFLKEKGGLDGLVWSPKRQEILDLWVYHTQGYFPDWQNYTP', init_coverage, pop)
    write_coverage_to_file("test_coverage_file.txt")


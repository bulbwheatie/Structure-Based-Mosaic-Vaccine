from pprint import pprint
import random

""" Important notes:

Functions to be called once at the very beginning of the program:
--calc_pop_epitope_freq()
--calc_single_freq() """

epitope_length = 9
possible_mutations = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
population_epitope_freq = {} # epitope string : frequency in population
population_num_epitopes = 0.0 # Total number of epitopes (non-distinct) in all the population sequences
population_top_single_freq = {}
hard_epitope_coverage = {} # epitope : coverage (according to Fisher)
        
def choose_mutation(mosaic_seq, init_coverage, population):
    """ Chooses a mutation by iterating though each position in the mosaic sequence and choosing, from the
    most frequently occuring 3 grams, what point mutations may increase coverage. """
    
    # Iterate through the rest of the amino acids in the middle
    # For each position, create some number of mutations, the most probable mutations conditioned
    # on the neighbors.  If the coverage increases for any of these mutations, then add them to a set
    # of possibilities to then choose from later.  This is a way to narrow the search space and potentially
    # get more coverage increases at each step (compared to random sampling)    
    top_choices = [] # (position, letter, coverage)
    rand_start = int(random.random() * epitope_length)
    for i in xrange(len(mosaic_seq)):
        # Look up neighbors.  If no neighbors, don't do anything (we have random sampling as a backup)
        max_mutations = 2
        num_mutations = min(max_mutations, len(population_top_single_freq[i]))
        mutation_choices = [population_top_single_freq[i][m][0] for m in xrange(num_mutations)]

        print mutation_choices
        for mutation_choice in mutation_choices:
            if mutation_choice != mosaic_seq[i]: # Only test if mutation is different than original sequence
                mutated_sequence = mosaic_seq[:i] + mutation_choice
                if i < len(mosaic_seq) - 1:
                    mutated_sequence += mosaic_seq[i+1:]
                curr_coverage = coverage(mutated_sequence)
                print curr_coverage
                if curr_coverage > init_coverage:
                    top_choices.append((i, mutation_choice, curr_coverage))

    if len(top_choices) == 0:
        return (-1, "-", init_coverage) # Flag that no mutations were found.
    else:
        # Introduce a degree of randomness here by choosing from the top 3 coverage-increasing mutations uniformly
        num_top_choices_considered = 3
        top_choices = sorted(top_choices, key=lambda x: x[2], reverse=True)
        num_considered = min(len(top_choices), num_top_choices_considered)
        return top_choices[int(random.random() * num_considered)]

def random_mutation(sequence):
    position = int(random.random() * len(sequence))
    amino_acid = possible_mutations[int(random.random() * 20)]
    return (position, amino_acid)

def num_epitopes_in_mosaic(mosaic, pop, eq_tolerance = 1, min_epitope_freq = 1):
    """ Returns the number of distinct population epitopes (with min_epitope_freq) that occur in the mosaic (with
        equality tolerance parameter) """
    epitopes_providing_coverage = set()
    for mosaic_epitope_start_i in xrange(len(mosaic) - epitope_length + 1):
        curr_mos_epi = mosaic[mosaic_epitope_start_i:mosaic_epitope_start_i + epitope_length]
        for key in population_epitope_freq:
            # Figure out if this particular population epitope is a near-match
            num_missed = 0
            for aa_i in xrange(epitope_length):
                if curr_mos_epi[aa_i] != key[aa_i]:
                    num_missed += 1
            if num_missed <= eq_tolerance and population_epitope_freq[key] >= min_epitope_freq:
                epitopes_providing_coverage.add(curr_mos_epi)
    return len(epitopes_providing_coverage)

def frac_pop_seq_covered(mosaic, pop, tolerance = False, coverage_thresh = 1):
    """ Returns the fraction of natural population sequences that have at least
        'coverage_thresh' epitopes covered by mosaic """
    # Construct set of mosaic epitopes to reduce runtime
    mosaic_epitopes = set()
    for mos_epi_start_i in xrange(len(mosaic) - epitope_length + 1):
        curr_epi = mosaic[mos_epi_start_i:mos_epi_start_i + epitope_length]
        mosaic_epitopes.add(curr_epi)
        if tolerance:
            # Deal with adding 9 * 19 more epitopes here
            pass

    # Look though each population sequence in turn
    num_pop_seq_covered = 0.0
    for pop_seq in pop:
        curr_match_count = 0
        for pop_epi_start_i in xrange(len(pop_seq) - epitope_length + 1):
            curr_pop_epi = pop_seq[pop_epi_start_i:pop_epi_start_i + epitope_length]
            if curr_pop_epi in mosaic_epitopes:
                curr_match_count += 1
                if curr_match_count >= coverage_thresh:
                    break
        if curr_match_count >= coverage_thresh:
            num_pop_seq_covered += 1.0
    return num_pop_seq_covered / len(pop)

def calc_single_freq(pop_seqs):
    """ Iterates through all aligned population sequences and calculates the amino acid frequence for each position"""
    global population_top_single_freq
    population_single_freq = {}
    for seq in pop_seqs:
        for i in xrange(len(seq)):
            amino_acid = seq[i]
            if i not in population_single_freq:
                population_single_freq[i] = {}
            if amino_acid not in population_single_freq[i]:
                population_single_freq[i][amino_acid] = 0.0
            population_single_freq[i][amino_acid] += 1.0

    for position in population_single_freq:
        position_probs = [(amino_acid, population_single_freq[position][amino_acid]) for amino_acid in population_single_freq[position]]
        population_top_single_freq[position] = sorted(position_probs, key=lambda x:x[1], reverse = True)

def calc_pop_epitope_freq(pop_seqs):
    """ Iterates through all population sequences and generates a dictionary of epitopes and frequency counts """
    global population_num_epitopes
    for seq in pop_seqs:
        for start_i in xrange(len(seq) - epitope_length + 1):
            curr_epitope = seq[start_i:start_i + epitope_length]
            if curr_epitope not in population_epitope_freq:
                population_epitope_freq[curr_epitope] = 0.0
            population_epitope_freq[curr_epitope] += 1.0
            
    # Sum frequencies of all epitopes
    for key in population_epitope_freq:
        population_num_epitopes += population_epitope_freq[key]

def coverage(mosaic_seq, threshold = 0.0):
    """ Iterate through mosaic epitopes.  For each key in the global population epitopes dictionary,
    add fractional coverage of that epitope * frequency.  Then at the end, divide by total number of epitopes.
    Range: 0 to 1 (for a specific epitope).  The coverage is the sum of sliding windows of epitopes in the mosaic,
    so it can and likely will be over 1 for longer sequences.
    """
    global population_num_epitopes
    total_coverage_score = 0.0
    for mosaic_start_i in xrange(len(mosaic_seq) - epitope_length + 1):
        for key in population_epitope_freq:
            if population_epitope_freq[key] >= threshold:
                epitope_coverage_score = 0.0
                for aa_i in xrange(epitope_length):
                    if mosaic_seq[mosaic_start_i + aa_i] == key[aa_i]:
                        epitope_coverage_score += 1.0
                epitope_coverage_score /= epitope_length
                total_coverage_score += epitope_coverage_score * population_epitope_freq[key]
    total_coverage_score /= population_num_epitopes
    return total_coverage_score

def fisher_coverage(mosaic_seq, population_seqs, threshold = 50):
    """ Returns coverage score for a mosaic sequence based on the population using Fisher's metric. """
    coverage = 0
    already_seen = set() # Avoids double counting coverage (within a mosaic protein)
    for start_i in xrange(len(mosaic_seq) - epitope_length + 1):
        curr_epitope = mosaic_seq[start_i:start_i + epitope_length]
        if curr_epitope in already_seen or curr_epitope not in population_epitope_freq:
            continue
        elif population_epitope_freq[curr_epitope] < threshold:
            continue
        else:
            already_seen.add(curr_epitope)

        if curr_epitope in hard_epitope_coverage:
            coverage += hard_epitope_coverage[curr_epitope]
        else:
            # Need to calculate coverage by iterating through all sequences.
            # In each iteration, calculate number of distinct epitopes.
            # If epitope in natural sequence, coverage = 1.0 / number of epitopes in sequence
            hard_epitope_coverage[curr_epitope] = 0.0
            for seq in population_seqs:
                curr_natural_seq_epitopes = set()
                epitope_match = 0.0
                for seq_start_i in xrange(len(seq) - epitope_length + 1):
                    if seq[seq_start_i:seq_start_i + epitope_length] == curr_epitope:
                        epitope_match = 1.0
                    curr_natural_seq_epitopes.add(seq[seq_start_i:seq_start_i + epitope_length])
                hard_epitope_coverage[curr_epitope] += epitope_match / len(curr_natural_seq_epitopes)

            hard_epitope_coverage[curr_epitope] /= len(population_seqs)
            coverage += hard_epitope_coverage[curr_epitope]
    return coverage

def read_fasta_file(fasta_file, start_i, end_i, aligned = True):
    sequences = []
    curr_seq = ""
    for line in open(fasta_file, 'r').readlines():
        if '>' in line:
            if len(curr_seq) != 0:
                sequences.append(curr_seq[start_i:end_i])
                curr_seq = ''
        else:
            curr_seq += line.replace('\n', '').strip()

    if not aligned:
        for i in xrange(len(sequences)):
            sequences[i] = sequences[i].replace('-', '')

    return sequences

# Tester function
if __name__ == "__main__":
    # Initialize data
    pop_gag = read_fasta_file('./data/HIV-1_gag.fasta', 343, 414, aligned=True)
    calc_pop_epitope_freq(pop_gag)
    calc_single_freq(pop_gag)
    print population_top_single_freq
    
    v1v2_seq = "SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"
    v1v2_start_i = 343
    v1v2_end_i = 414
    init_coverage = coverage(v1v2_seq)
    print init_coverage
    print choose_mutation(v1v2_seq, init_coverage, pop_gag)
    print fisher_coverage(v1v2_seq, pop_gag)
    print "num epitopes in mosaic: ", num_epitopes_in_mosaic(v1v2_seq, pop_gag, eq_tolerance = 1, min_epitope_freq = 1)
    print "frac pop covered: ", frac_pop_seq_covered(v1v2_seq, pop_gag, tolerance = False, coverage_thresh = 20)

    nef_seq = 'AWLEAQEEEEVGFPVTPQVPLRPMTYKAAVDLSHFLKEKGGLEGLIHSQRRQDILDLWIYHTQGYFPDWQNYTPGPGIRYPLTFGWCYKLVPVEPEKLEEANKDDPEREVLEWRFDSRLAFHHMARELHPEYFKNA'
    
    #output_seq = 'SILKEGPGPAEPLRQILGLGLEEEEEEEAEEEEEEEELEALLVKGAPGLSKTLLKALGEGATLEEALGGHQ'
    #print "output mosaic num epitopes: ", num_epitopes_in_mosaic(output_seq, pop, eq_tolerance = 1, min_epitope_freq = 1)
    #print "frac pop covered: ", frac_pop_seq_covered(output_seq, pop, tolerance = False, coverage_thresh = 20)

################################## DEPRECATED FUNCTIONS ####################################

def d_write_coverage_to_file(filename, cov_type = 'soft'):
    with open(filename, 'w') as f:
        dict_to_write = None
        if cov_type == 'soft':
            dict_to_write = soft_epitope_coverage
        else:
            dict_to_write = hard_epitope_coverage
        for epitope in dict_to_write:
            f.write(epitope + ',' + str(dict_to_write[epitope]) + '\n')
    
def d_read_coverage_from_file(filename, cov_type = 'soft'):
    lines = open(filename, 'r').readlines()
    for l in lines:
        epitope, coverage = l.split(',')
        if cov_type == 'soft':
            soft_epitope_coverage[epitope] = float(coverage)
        else:
            hard_epitope_coverage[epitope] = float(coverage)

def d_get_location_probabilities(mosaic_seq):
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

def d_coverage(mosaic_seq, population_seqs, t="soft"):
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
                sequence_match_overall = 0.0
                for seq_start_i in xrange(len(seq) - epitope_length + 1):
                    epitope_match = 0.0
                    if t == 'hard':
                        if seq[seq_start_i:seq_start_i + epitope_length] == curr_epitope:
                            epitope_match = 1.0
                            break
                    else: # type == 'soft'
                        aa_position_match_count = 0.0
                        for aa_position in xrange(epitope_length):
                            if curr_epitope[aa_position] == seq[seq_start_i + aa_position]:
                                aa_position_match_count += 1.0
                        perfect_match_weight = 1.0
                        if aa_position_match_count == epitope_length:
                            epitope_match += perfect_match_weight # To weight perfect epitope matches highery
                        else:
                            epitope_match += aa_position_match_count / epitope_length

                    sequence_match_overall += epitope_match
                    curr_natural_seq_epitopes.add(seq[seq_start_i:seq_start_i + epitope_length])
                if t == 'soft':
                    soft_epitope_coverage[curr_epitope] += sequence_match_overall / len(seq)
                else:
                    hard_epitope_coverage[curr_epitope] += sequence_match_overall / len(curr_natural_seq_epitopes)
                print sequence_match_overall, len(curr_natural_seq_epitopes), len(seq)

            print "coverage, then soft_epitope_coverage[]", coverage, soft_epitope_coverage[curr_epitope]
            if t == 'soft':
                soft_epitope_coverage[curr_epitope] /= len(population_seqs)
                coverage += soft_epitope_coverage[curr_epitope]
            else:
                hard_epitope_coverage[curr_epitope] /= len(population_seqs)
                coverage += hard_epitope_coverage[curr_epitope] 
            
    return coverage

def d_get_all_sequences(aligned = True):
    """ This function returns a list of all fasta sequences for HIV from
        the global list, fasta_filenames.

        Remove '-' placeholders that code for empty spaces in the FASTA files """
    
    all_seqs = []
    for file in fasta_filenames:
        all_seqs.extend(read_fasta_file(file, aligned))
    return all_seqs

def d_calc_aa_ngrams(pop_seqs):
    """ Calculate the conditional probability for each amino acid given it's left and right neighbors.
    This is an effort to narrow the search space for mutations to more common sequence, and hopefully
    more common epitopes.  """
    
    aa_ngrams = {} # prior & after : middle : count
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

    # Divide frequencies by counts to get conditional probabilities
    for neighbor in aa_ngrams:
        curr_sum = 0.0
        for middle in aa_ngrams[neighbor]:
            curr_sum += aa_ngrams[neighbor][middle]

        for middle in aa_ngrams[neighbor]:
            aa_ngrams[neighbor][middle] /= curr_sum

    # Sort the ngrams based on probability for future use
    for neighbor in aa_ngrams:
        middle_probabilities = [(middle, aa_ngrams[neighbor][middle]) for middle in aa_ngrams[neighbor]]
        aa_top_ngrams[neighbor] = sorted(middle_probabilities, key=lambda x: x[1], reverse=True)[:5]

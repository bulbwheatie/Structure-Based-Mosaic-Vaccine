import random, math

""" Important notes:

Functions to be called once at the very beginning of the program:
--calc_pop_epitope_freq()
--calc_single_freq() """

epitope_length = 9
possible_mutations = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
population_epitope_freq_aligned = {} # mosaic position index: ALIGNED epitope string : frequency in population
population_epitope_freq_unaligned = {} # UNALIGNED epitope string : frequency in population
population_num_epitopes = 0.0 # Total number of epitopes (non-distinct) in all the population sequences
population_top_single_freq = {} # ALIGNED position : (amino acid, frequency) Stores the amino acid frequences at each position
population_top_nmer_freq = {} # ALIGNED n(length of nmer) : position : amino acid : frequency (dict of dicts of nmer frequences)
hard_epitope_coverage = {} # UNALIGNED epitope : coverage (according to Fisher)
soft_epitope_coverage = {} # mosaic position index: UNALIGNED epitope : coverage (according to us usign a particular weight function)
num_population_sequences = None


# Functions that weight partial matches (for our own definition of convergence)
squared_denominator = float(epitope_length ** 2)
squared_numerator = [float(i ** 2) for i in xrange(epitope_length + 1)]
def squared_weight(num_matches):
    return squared_numerator[num_matches] / squared_denominator

exponential_denominator = float(2 ** epitope_length - 1)
exponential_numerator = [float(2 ** i - 1) for i in xrange(epitope_length + 1)]
def exponential_weight(num_matches):
    return exponential_numerator[num_matches] / exponential_denominator

def choose_point_mutation(mosaic_seq, init_coverage, max_mutations_per_position = 2, position_mutation_probability = 0.5, allow_insertions_deletions = False, weight_func = exponential_weight, coverage_temperature = .005):
    """ Chooses a mutation by iterating though each position in the mosaic sequence and choosing, from the
    most frequently occuring 3 grams, what point mutations may increase coverage. """
    
    global num_population_sequences
    top_choices = [] # (position, letter, coverage)

    # For each position, create some number of mutations, the most probable mutations conditioned
    # on the neighbors.  If the coverage increases for any of these mutations, then add them to a set
    # of possibilities to then choose from later.  This is a way to narrow the search space and potentially
    # get more coverage increases at each step (compared to random sampling)
    for i in xrange(len(mosaic_seq)):
        mutation_choices = []
        # Only look at mutations for around half the positions to save time
        if random.random() < position_mutation_probability:
            # Look up neighbors.  If no neighbors, don't do anything (we have random sampling as a backup)
            num_mutations = min(max_mutations_per_position, len(population_top_single_freq[i]))

            # Choose 'num_mutations' mutations proportional to their frequency in the population
            # Generate random integers between 0 and sum(freq), then use the mutation choice that
            # corresponds to that number.  To avoid repeats, remove previous mutation frequency
            # from the random number range
            rand_num_range = num_population_sequences
            for m in xrange(num_mutations):
                leftover = int(random.random() * rand_num_range)
                count = 0
                last_mutation = None
                last_mutation_freq = None
                while leftover >= 0:
                    if population_top_single_freq[i][count][0] not in mutation_choices:
                        leftover -= population_top_single_freq[i][count][1]
                        last_mutation = population_top_single_freq[i][count][0]
                        last_mutation_count = population_top_single_freq[i][count][1]
                    count += 1
                if (allow_insertions_deletions or
                    (last_mutation == "-" and mosaic_seq[i] == "-") or (last_mutation != "-" and mosaic_seq[i] != "-")):
                    mutation_choices.append(last_mutation)
                rand_num_range -= last_mutation_count

        print "Considering the following mutation choice(s) at position " + str(i) + ":", mutation_choices
        for mutation_choice in mutation_choices:
            if mutation_choice != mosaic_seq[i]: # Only test if mutation is different than original sequence
                mutated_sequence = update_seq_string(mosaic_seq, mutation_choice, i)
                curr_coverage = coverage(mutated_sequence, weight_func = weight_func)
                print "The coverage (if making mutation" + mutation_choice + "at position " + str(i) + ") would be ", curr_coverage
                improvement_threshold = .005
                if curr_coverage >= init_coverage + improvement_threshold:
                    top_choices.append((i, mutation_choice, curr_coverage))
                else:
                    prob_accept = math.pow(math.e, (curr_coverage - init_coverage - improvement_threshold) / coverage_temperature)
                    if random.random() < prob_accept:
                        top_choices.append((i, mutation_choice, curr_coverage))

    if len(top_choices) == 0:
        return (-1, "-", init_coverage) # Flag that no mutations were found.
    else:
        # Introduce a degree of randomness here by choosing from the top 3 coverage-increasing mutations uniformly
        num_top_choices_considered = 3
        top_choices = sorted(top_choices, key=lambda x: x[2], reverse=True)
        #num_considered = min(len(top_choices), num_top_choices_considered)
        num_considered = len(top_choices)
        return top_choices[int(random.random() * num_considered)]

def choose_n_sub_mutation(mosaic_seq, init_coverage, pop, mut_length = 2, max_mutations_per_position = 2, position_mutation_probability = 0.5, weight_func = exponential_weight, coverage_temperature = 0.005):
    """ Choose a substution mutation of length 'mut_length' that represents no insertion/deletions.
        Chooses the substitution based on the highest_frequency nmer starting at each position, ruling
        out substitutions that:
        (1) have a gap in the mosaic sequence where the selected substitution doesn't have a gap.
        (2) have an amino acid where the selected subtitution has a gap.

        This function returns ('mut_length' - the number of matching amino acids) mutations """
    global population_top_nmer_freq # n(length of nmer) : position : amino acids : frequency (dict of dicts of nmer frequences)
    # We need to create a dictionary of mutations of length 'max_mutations_per_position'.  This is used to choose
    # more frequently occuring chunks to place in the mosaic.
    # Compiling a dicitonary only needs to be done once, so check if mut_length in the dictionary 'population_top_nmer_freq'
    if mut_length not in population_top_nmer_freq:
        # If not, we need to populate the dictionary.
        # Iterate through population sequence positions to calculate frequencies of each nmer at each position.
        population_top_nmer_freq[mut_length] = {}
        for pop_seq in pop:
            for aa_i in xrange(len(pop_seq) - mut_length + 1):
                if aa_i not in population_top_nmer_freq[mut_length]:
                    population_top_nmer_freq[mut_length][aa_i] = {}
                curr_nmer = pop_seq[aa_i:aa_i + mut_length]
                if curr_nmer not in population_top_nmer_freq[mut_length][aa_i]:
                    population_top_nmer_freq[mut_length][aa_i][curr_nmer] = 0.0
                population_top_nmer_freq[mut_length][aa_i][curr_nmer] += 1.0

        # Collapse the amino acids : frequency dictionary into a tuple to facilitate sorting by frequency
        for aa_i in population_top_nmer_freq[mut_length]:
            temp_list = []
            for nmer in population_top_nmer_freq[mut_length][aa_i]:
                temp_list.append((nmer, population_top_nmer_freq[mut_length][aa_i][nmer]))
            population_top_nmer_freq[mut_length][aa_i] = sorted(temp_list, key=lambda x: x[1], reverse=True)

    # Follow a similar path to the single point mutation chooser, though with more checks against
    # substitution/insertion/deletion mixes.
    top_choices = [] # (position, letters, coverage)
    for i in xrange(len(mosaic_seq) - mut_length + 1):
        mutation_choices_unpruned = []
        if random.random() < position_mutation_probability:
            num_mutations = min(max_mutations_per_position, len(population_top_nmer_freq[mut_length][i]))

            # Choose 'num_mutations' mutations proportional to their frequency in the population
            rand_num_range = num_population_sequences
            for m in xrange(num_mutations):
                leftover = int(random.random() * rand_num_range)
                count = 0
                last_mutation = None
                last_mutation_freq = None
                while leftover >= 0:
                    if population_top_nmer_freq[mut_length][i][count][0] not in mutation_choices_unpruned:
                        leftover -= population_top_nmer_freq[mut_length][i][count][1]
                        last_mutation = population_top_nmer_freq[mut_length][i][count][0]
                        last_mutation_count = population_top_nmer_freq[mut_length][i][count][1]
                    count += 1
                mutation_choices_unpruned.append(last_mutation)
                rand_num_range -= last_mutation_count

        mutation_choices = []
        # Iterate through mutation choices, removing any that result in insertions/deletions
        for mutation_choice in mutation_choices_unpruned:
            valid = True
            for aa_i in xrange(len(mutation_choice)):
                if ((mutation_choice[aa_i] == "-" and mosaic_seq[i + aa_i] != "-")
                    or (mutation_choice[aa_i] != "-" and mosaic_seq[i + aa_i] == "-")):
                    valid = False
            if valid:
                mutation_choices.append(mutation_choice)
        print "Considering the following mutation choice(s) at position " + str(i) + ":", mutation_choices
        
        for mutation_choice in mutation_choices:
            if mutation_choice != mosaic_seq[i:i + mut_length]: # Only test if mutation is different than original sequence
                mutated_sequence = update_seq_string(mosaic_seq, mutation_choice, i)
                curr_coverage = coverage(mutated_sequence, weight_func = weight_func)
                improvement_threshold = 0.005
                print curr_coverage
                if curr_coverage > init_coverage + improvement_threshold:
                    top_choices.append((i, mutation_choice, curr_coverage))
                else:
                    prob_accept = math.pow(math.e, (curr_coverage - init_coverage - improvement_threshold) / coverage_temperature)
                    if random.random() < prob_accept:
                        top_choices.append((i, mutation_choice, curr_coverage))


    if len(top_choices) == 0:
        return [(-1, "-", init_coverage)] # Flag that no mutations were found.
    else:
        # Introduce a degree of randomness here by choosing from the top 3 coverage-increasing mutations uniformly
        num_top_choices_considered = 3
        top_choices = sorted(top_choices, key=lambda x: x[2], reverse=True)
        #num_considered = min(len(top_choices), num_top_choices_considered)
        num_considered = len(top_choices)
        final_choice = top_choices[int(random.random() * num_considered)]

        # Reformat final_choice as list of mutations
        formatted_mutation = []
        for aa_i in xrange(mut_length):
            if final_choice[1][aa_i] != mosaic_seq[final_choice[0] + aa_i]:
                formatted_mutation.append((final_choice[0] + aa_i, final_choice[1][aa_i], final_choice[2]))
        return formatted_mutation

def calc_single_freq(pop_seqs):
    """ Iterates through all aligned population sequences and calculates the amino acid frequence for each position"""
    global population_top_single_freq
    global num_population_sequences
    num_population_sequences = len(pop_seqs)
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
    """ Iterates through all population sequences and generates a dictionary of epitopes and frequency counts
    PASS IN ALIGNED POPULATION """
    global population_epitope_freq_aligned
    global population_epitope_freq_unaligned
    global population_num_epitopes

    for temp_i in xrange(len(pop_seqs[0]) - epitope_length + 1):
        population_epitope_freq_aligned[temp_i] = {}
    
    for seq in pop_seqs:
        # Aligned
        for start_i in xrange(len(seq) - epitope_length + 1):
            curr_epitope = seq[start_i:start_i + epitope_length]
            if curr_epitope not in population_epitope_freq_aligned[start_i]:
                population_epitope_freq_aligned[start_i][curr_epitope] = 0.0
            population_epitope_freq_aligned[start_i][curr_epitope] += 1.0
        # Unaligned
        unaligned_seq = seq.replace("-", "")
        for start_i in xrange(len(unaligned_seq) - epitope_length + 1):
            curr_epitope = unaligned_seq[start_i:start_i + epitope_length]
            if curr_epitope not in population_epitope_freq_unaligned:
                population_epitope_freq_unaligned[curr_epitope] = 0.0
            population_epitope_freq_unaligned[curr_epitope] += 1.0
            
    # Sum frequencies of all epitopes
    for key in population_epitope_freq_unaligned:
        population_num_epitopes += population_epitope_freq_unaligned[key]

def coverage(mosaic_seq, threshold = 0.0, weight_func = exponential_weight):
    """ Iterate through mosaic epitopes.  For each key in the global population epitopes dictionary,
    add fractional coverage of that epitope * frequency.  Then at the end, divide by total number of epitopes.
    Range: 0 to 1 (for a specific epitope).  The coverage is the sum of sliding windows of epitopes in the mosaic,
    so it can and likely will be over 1 for longer sequences.
    """
    global population_num_epitopes
    global soft_epitope_coverage

    if len(soft_epitope_coverage) == 0:
        for temp_i in xrange(len(mosaic_seq) - epitope_length + 1):
            soft_epitope_coverage[temp_i] = {}
            
    #mosaic_seq = mosaic_seq.replace("-", "")
    total_coverage_score = 0.0
    for mosaic_start_i in xrange(len(mosaic_seq) - epitope_length + 1):
        curr_mosaic_epi = mosaic_seq[mosaic_start_i:mosaic_start_i + epitope_length]
        if '-' in curr_mosaic_epi:
            continue
        if curr_mosaic_epi in soft_epitope_coverage[mosaic_start_i]:
            # We've cached the result; let's pull it out to reduce computation time
            total_coverage_score += soft_epitope_coverage[mosaic_start_i][curr_mosaic_epi]
        else:
            # First time we're seeing this epitope, so let's calculate and then cache
            for key in population_epitope_freq_aligned[mosaic_start_i]:
                if population_epitope_freq_aligned[mosaic_start_i][key] >= threshold:
                    epitope_coverage_score = 0
                    for aa_i in xrange(epitope_length):
                        if mosaic_seq[mosaic_start_i + aa_i] == key[aa_i]:
                            epitope_coverage_score += 1
                    epitope_coverage_score = weight_func(epitope_coverage_score)
                    if curr_mosaic_epi not in soft_epitope_coverage[mosaic_start_i]:
                        soft_epitope_coverage[mosaic_start_i][curr_mosaic_epi] = 0.0
                    soft_epitope_coverage[mosaic_start_i][curr_mosaic_epi] += epitope_coverage_score * population_epitope_freq_aligned[mosaic_start_i][key]
            total_coverage_score += soft_epitope_coverage[mosaic_start_i][curr_mosaic_epi]
    total_coverage_score /= population_num_epitopes
    return total_coverage_score

def fisher_coverage(mosaic_seq, population_seqs, threshold = 50):
    """ Returns coverage score for a mosaic sequence based on the population using Fisher's metric.
        Note: Population sequences should be ALIGNED! """
    coverage = 0
    mosaic_seq = mosaic_seq.replace("-", "")
    already_seen = set() # Avoids double counting coverage (within a mosaic protein)
    for start_i in xrange(len(mosaic_seq) - epitope_length + 1):
        curr_epitope = mosaic_seq[start_i:start_i + epitope_length]
        if curr_epitope in already_seen or curr_epitope not in population_epitope_freq_unaligned:
            continue
        elif population_epitope_freq_unaligned[curr_epitope] < threshold:
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
                seq_ungapped = seq.replace("-", "")
                curr_natural_seq_epitopes = set()
                epitope_match = 0.0
                for seq_start_i in xrange(len(seq_ungapped) - epitope_length + 1):
                    if seq_ungapped[seq_start_i:seq_start_i + epitope_length] == curr_epitope:
                        epitope_match = 1.0
                    curr_natural_seq_epitopes.add(seq_ungapped[seq_start_i:seq_start_i + epitope_length])
                hard_epitope_coverage[curr_epitope] += epitope_match / len(curr_natural_seq_epitopes)

            hard_epitope_coverage[curr_epitope] /= len(population_seqs)
            coverage += hard_epitope_coverage[curr_epitope]
    return coverage

def update_seq_string(mosaic_seq, mutation_choice, pos):
    """ Simple utility function to update and return a mosaic after a point mutation.
    pos is the popition at which the 'mutation_choice' is being substituted into the sequence (mosaic_seq) """
    
    mutated_sequence = mosaic_seq[:pos] + mutation_choice
    if pos < len(mosaic_seq) - len(mutation_choice):
        mutated_sequence += mosaic_seq[pos + len(mutation_choice):]
    return mutated_sequence

def prune_population_seqs(mosaic_seq, pop_seqs):
    """ Removes all gaps in the population sequence where the mosaic sequence has gaps. """
    i_to_remove = []
    for i in xrange(len(mosaic_seq) - 1, -1, -1):
        if mosaic_seq[i] == '-':
            i_to_remove.append(i)

    for pop_i in xrange(len(pop_seqs)):
        for ii in i_to_remove:
            if ii != len(pop_seqs[pop_i]) - 1:
                pop_seqs[pop_i] = pop_seqs[pop_i][:ii] + pop_seqs[pop_i][(ii + 1):]
            else:
                pop_seqs[pop_i] = pop_seqs[pop_i][:ii]
                
def read_fasta_file(fasta_file, start_i, end_i, aligned = True):
    """ Reads the population sequences that come in FASTA file format.
    start_i and end_i denote the starting and ending indices for the mosaic beginning/end points with respect
    to the populations equences. """
    
    sequences_all = []
    curr_seq = ""
    for line in open(fasta_file, 'r').readlines():
        if '>' in line:
            if len(curr_seq) != 0:
                sequences_all.append(curr_seq[start_i:end_i])
                curr_seq = ''
        else:
            curr_seq += line.replace('\n', '').strip()

    # Remove sequences that aren't the same length.  There are some oddball nef ones that aren't the same length
    sequences = []
    for seq in sequences_all:
        if len(seq) == end_i - start_i:
            sequences.append(seq)

    if not aligned:
        for i in xrange(len(sequences)):
            sequences[i] = sequences[i].replace('-', '')

    return sequences

def find_ins_del_positions(seqs):
    gap_counts = [0.0] * len(seqs[0])
    for seq in seqs:
        for aa_i in xrange(len(seq)):
            if seq[aa_i] == '-':
                gap_counts[aa_i] += 1.0

    possible_mutation_sites = [i for i in xrange(len(seqs[0])) if (gap_counts[i] > 200 and gap_counts[i] < len(seqs) - 200)]
    possible_gap_counts = [gap_counts[i] for i in possible_mutation_sites]
    print zip(possible_mutation_sites, possible_gap_counts)
    print len(seqs), gap_counts
    return possible_mutation_sites














################# Functions that aren't currently used but may be useful for later #####################
# Denoted with the word 'deprecated' to avoid confusion for outside readers trying to understand how the program works.

def deprecated_random_mutation(sequence):
    """ Chooses a random point mutation in the sequence.  Deprecated, but may be useful for the future. """
    position = int(random.random() * len(sequence))

    while (sequence[position] == "-"):
        position = int(random.random() * len(sequence))

    amino_acid = possible_mutations[int(random.random() * 20)]

    #Calculate the position in the pose
    pose_position = 0
    for i in xrange(position + 1):
        if (sequence[i] != "-"):
            pose_position +=1

    #Pose position is zero indexed
    return (position, pose_position + 1, amino_acid)

def deprecated_num_epitopes_in_mosaic(mosaic, pop, eq_tolerance = 1, min_epitope_freq = 1):
    """ Returns the number of distinct population epitopes (with min_epitope_freq) that occur in the mosaic (with
        equality tolerance parameter) """
    mosaic = mosaic.replace("-", "")
    epitopes_providing_coverage = set()
    for mosaic_epitope_start_i in xrange(len(mosaic) - epitope_length + 1):
        curr_mos_epi = mosaic[mosaic_epitope_start_i:mosaic_epitope_start_i + epitope_length]
        for key in population_epitope_freq_unaligned:
            # Figure out if this particular population epitope is a near-match
            num_missed = 0
            for aa_i in xrange(epitope_length):
                if curr_mos_epi[aa_i] != key[aa_i]:
                    num_missed += 1
            if num_missed <= eq_tolerance and population_epitope_freq_unaligned[key] >= min_epitope_freq:
                epitopes_providing_coverage.add(curr_mos_epi)
    return len(epitopes_providing_coverage)

def deprecated_frac_pop_seq_covered(mosaic, pop, tolerance = False, coverage_thresh = 1):
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


# Tester function
if __name__ == "__main__":
    # Find sites for ins/del
    v1v2_seq = 'VKLTPLCVTLQCTNVTNNITD-------------------------------------DMRGELKN----CSFNM-T-TE--LRD-KK-QKV-YSLF-YRLDVVQINENQGNRSNNS------------------------------------------NKEYRLI---NCNTSAI-T---QA'
    v1v2_start_i = 171
    v1v2_end_i = 354
    pop_env = read_fasta_file('./data/HIV-1_env.fasta', v1v2_start_i, v1v2_end_i, aligned=True)
    insdel_sites = find_ins_del_positions(pop_env)
    for site in insdel_sites:
        print site, v1v2_seq[site]

    nef_start_i = 111
    nef_end_i = 357
    pop_nef_aligned = read_fasta_file('./data/HIV-1_nef.fasta', nef_start_i, nef_end_i, aligned = True)
    #find_ins_del_positions(pop_nef_aligned)

    print "v1v2 loop diagnostics"
    consensus = 'VKLTPLCVTLNCTDVNNTNTTMEKGEIKNCSFNITTEIRDKVQKEYALFYKLDVVPIDNNNNSNNSSSNTSYRLINCNTSVITQA'
    prune_population_seqs(v1v2_seq, pop_env)
    calc_pop_epitope_freq(pop_env)
    calc_single_freq(pop_env)

    v1v2_seq = v1v2_seq.replace("-","")
    with open('temp_popv1v2.fasta', 'w') as f:
        for seq in pop_env:
            f.write(seq + '\n')

    init_coverage = coverage(v1v2_seq, weight_func = squared_weight)
    consensus_coverage = coverage(consensus, weight_func = squared_weight)
    init_hard = fisher_coverage(v1v2_seq, pop_env)
    consensus_hard = fisher_coverage(consensus, pop_env)

    print "Init coverage: ", init_coverage, consensus_coverage, init_hard, consensus_hard
    #print choose_point_mutation(v1v2_seq, init_coverage, allow_insertions_deletions = True)
    #choose_n_sub_mutation(v1v2_seq, init_coverage, pop_env, mut_length = 2, max_mutations_per_position = 1)

    """print "gag loop diagnostics"
    gag_start_i = 343
    gag_end_i = 414
    pop_gag = read_fasta_file('./data/HIV-1_gag.fasta', gag_start_i, gag_end_i, aligned=True)
    calc_pop_epitope_freq(pop_gag)
    calc_single_freq(pop_gag)
    gag_seq = "SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"
    init_coverage = coverage(gag_seq, weight_func=exponential_weight)
    print population_top_nmer_freq
    #for i in xrange(11):
    #    print choose_n_sub_mutation(gag_seq, init_coverage, pop_gag, mut_length = 2, max_mutations_per_position = 3)
    print init_coverage
    print choose_point_mutation(gag_seq, init_coverage)
    print fisher_coverage(gag_seq, pop_gag)"""

    print
    print "Nef sequence diagnostics"
    #nef_seq = 'AWLEAQEEEEVGFPVTPQVPLRPMTYKAAVDLSHFLKEKGGLEGLIHSQRRQDILDLWIYHTQGYFPDWQNYTPGPGIRYPLTFGWCYKLVPVEPEKLEEANKDDPEREVLEWRFDSRLAFHHMARELHPEYFKNA'
    nef_seq = 'AWL--EA-QE-----E---E--E--VGFPVTPQVPLRPMTYKAAVDLSHFLKEKGGLEGLIHSQRRQDILDLWIYHTQGYFPDWQNYTPGPGIRYP-----------------------------------------------------------------LTFGWCYKLVPVEPEKLE-EANK---------------------------DDP-EREVLEWRFDSRLAFHHMARELHPEYF-KNA'
    consensus = 'AWLEAQEEEEVGFPVRPQVPLRPMTYKGAFDLSHFLKEKGGLEGLIYSQKRQDILDLWVYHTQGYFPDWQNYTPGPGIRYPLTFGWCFKLVPVDPEEVEEANEDDPEREVLMWKFDSRLAFRHMARELHPEYYKLK'

    prune_population_seqs(nef_seq, pop_nef_aligned)
    nef_seq = nef_seq.replace("-","")
    with open('popnef.fasta', 'w') as f:
        for seq in pop_nef_aligned:
            f.write(seq + '\n')
    pop_nef_unaligned = read_fasta_file('./data/HIV-1_nef.fasta', nef_start_i, nef_end_i, aligned = False)
    calc_pop_epitope_freq(pop_nef_aligned)
    calc_single_freq(pop_nef_aligned)
    print coverage(nef_seq, weight_func=squared_weight)    
    
    init_coverage = coverage(nef_seq, weight_func = squared_weight)
    consensus_coverage = coverage(consensus, weight_func = squared_weight)
    init_hard = fisher_coverage(nef_seq, pop_nef_aligned)
    consensus_hard = fisher_coverage(consensus, pop_nef_aligned)
    print init_coverage, consensus_coverage, init_hard, consensus_hard

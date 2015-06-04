"""
struct_utils.py 
------------------
MUST CALL rosetta.init()
MUST CALL initialize_struct_utils()
MUST CALL set_fLog()
...prior to calling functions in this file 
------------------

This file contains functions for performing mutations and actions related to 
structure. It maintains the connection between structure and aligned sequences
and translates output from utils.py into structural changes. 
"""

from rosetta import *
from utils import *
from optimizeStructure import *
import math
import sys

possible_mutations = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# Call zero_pose() to create native structure aligned for deleteion and insertion purposes


#Keep a current [MASTER] 7aligned sequence (contains gaps) 
sequence_archive = [0] 
pose_archive = [0] 
cover_archive = [0]#Initialize to non convergent values
prev_top_ten_cover = [0] * 10
same_top_ten_count = 0
archive_idx = 0
convergence_flag = False
intermediate_struct_counter = 0
name_space = ""
log = None
reject_count = 0
number_of_resets = 0
mutation_temp = 0.005
preset_idx = 0
preset_attempt = 1
consensus_sequence = ""

point_to_chunk_prob = 0.5

#When we make a mutation, iterate determine find correlating position in structure
def calculate_mutation_for_pose(sequence_master, position, amino_acid, count):
	""" Compute the corresponding indices from the master sequence for the pose object
	Returns the pose position for where to mutate, delete or insert before
	Returns the type of mutation to make
	If the mutation is invalid (count == -1), make a random mutation
	"""
	#Count the number of existing residues
	pose_position = 0
	curr_position = 0
	mutation_type = "point"

	#If the mutation chooser returned nothing, then choose a random one
	if (count < 0):
		(position, amino_acid) = random_mutation(sequence_master)

	while curr_position < position:
		#Residue exists in the structure
		if (sequence_master[curr_position] != "-"):
			pose_position += 1
		curr_position += 1

	if (sequence_master[position] == "-"):
		mutation_type = "insert"
	elif (amino_acid == "-"):
		mutation_type = "delete"

	return (pose_position + 1, mutation_type)

def make_mutation(pop, coverage_weight, mutation_length = 2, follow_consensus = False):
	"""
		(1) Gets the current sequence
		(2) Randomly chooses either a point or chunk mutation 
		(3) Attempts to make that mutation 
			(3.1) Success = break
			(3.2) Fail = attempt mutation of other type
				(3.2.1) Success = break
				(3.2.2) Fail = random til successs


	"""
	(tmp_pose, sequence, cover) = get_current_structure()

	
	mutation_type = ""
	if (tmp_pose == 0 or sequence == 0):
		#Rejected too many structures, terminate
		return (-1, -1, -1, -1, -1)

	pose = Pose()
	pose.assign(tmp_pose)

	if (follow_consensus):
		#TODO: Make a mutation according to the consensus sequence
		(pose, position, mutation) = mutate_preset_sequence(pose)
		return (pose, pose.sequence(), position, mutation, "PRESET")


	if (random.random() > point_to_chunk_prob):
		#Make a chunk mutation first
		(tmp_sequence, positions, mutations) = make_chunk_mutation(pose, sequence, cover, pop, coverage_weight, mutation_length) 
		mutation_type = "CHUNK"
		if (tmp_sequence == -1):
			#Make point mutation
			(tmp_sequence, positions, mutations) = make_point_mutation(pose, sequence, cover, coverage_weight) 
			mutation_type = "POINT"
			if (tmp_sequence == -1):
				#Don't Make a random mutation
				#(sequence, positions, mutations) = make_random_mutation(pose, sequence)
				mutation_type = "NONE"
				positions = "-"
				mutations = "-"
			else:
				sequence = tmp_sequence
		else:
			sequence = tmp_sequence
	else:
		#Make a point mutation first
		(tmp_sequence, positions, mutations) = make_point_mutation(pose, sequence, cover, coverage_weight) 
		mutation_type = "POINT"
		if (tmp_sequence == -1):
			#Make point mutation
			(tmp_sequence, positions, mutations) = make_chunk_mutation(pose, sequence, cover, pop, coverage_weight, mutation_length) 
			mutation_type = "CHUNK"
			if (tmp_sequence == -1):
				#Don't Make a random mutation
				#(sequence, positions, mutations) = make_random_mutation(pose, sequence)
				mutation_type = "NONE"
				positions = "-"
				mutations = "-"
			else:
				sequence = tmp_sequence
		else:
			sequence = tmp_sequence

	#Now sequence is our mutated sequence
	#Pose is our mutated pose

	return (pose, sequence, positions, mutations, mutation_type)

def make_point_mutation(pose, sequence, cover, coverage_weight):
	"""
	Make a point mutation by calling utils and then using the PyRosetta method
	If the mutation fails for any reason, return an error code of (-1, -1, -1)
	"""
	(position, mutation, cover) = choose_point_mutation(sequence, cover, allow_insertions_deletions = False, weight_func = coverage_weight, coverage_temperature=mutation_temp)
	log.write("POINT MUTATION: " + str(position) + " to " + mutation + "\n")
	if (position == -1):
		return (-1, -1, -1)
	(pose_position, mutation_type) = calculate_mutation_for_pose(sequence, position, mutation, 0)
	print mutation + " at " + str(position) + "\n"
	if ("-" in mutation and position > -1):
		print "ERROR: - in mutation choice\n"
		log.write("ERROR: - in mutation choice\n")
		return (-1, -1, -1) 
	#Check the type and make the appropriate point mutation
	if (mutation_type == "point" ):
		print "Point\n"
		try:
			mutate_residue(pose, pose_position, mutation)
		except:
			log.write("ERROR: Could not mutate structure \n")
			log.close()
			raise
		if (sequence[position] == pose.sequence()[pose_position - 1]):
			return (-1, -1, -1)
	elif (mutation_type == "insert"):
		dump_pdb(pose, "tmpRemodelPdb.pdb")
		(new_pdb, mutated_sequence) = insert_residue("tmpRemodelPdb.pdb", position, sequence)
		pose_from_pdb(pose, new_pdb)
		print "Insert\n"
		return (mutated_sequence, position, pose.sequence()[pose_position + 1])
	elif (mutation_type == "delete"):
		dump_pdb(pose, new_pdb)
		(new_pdb, mutated_sequence) = delete_residue("tmpRemodelPdb.pdb", position, sequence)
		print "Delete\n"
		return (mutated_sequence, position, "-")

	mutated_sequence = update_seq_string(sequence, mutation, position)

	log.write("POINT MUTATION: sequence =" + mutated_sequence + "\n")
	return (mutated_sequence, position, mutation)

def make_chunk_mutation(pose, sequence, cover, pop, coverage_weight, mutation_length):
	"""
	Make a chunk mutation by calling utils and then making
	each mutation as a point mutation. If any of the mutations fail
	revert the entire structure to before the chunk mutation
	"""
	mutations = choose_n_sub_mutation(sequence, cover, pop, mut_length = mutation_length, weight_func = coverage_weight, coverage_temperature=mutation_temp)
	log.write("CHUNK MUTATION: " + str(mutations)  + "\n")
	log_positions = []
	log_mutations = []
	print str(mutations) + "\n"
	# Keep a temporary pose in case not everything is successful
	test_pose = Pose()
	test_pose.assign(pose)

	mutated_sequence = sequence
	if (mutations is not None and mutations[0] != -1 and mutations[0][0] != -1):
		i = 0
		while i < len(mutations):
			if ("-" in mutations[i] and positions[i] > -1):
				print "ERROR: - in mutation choice\n"
				log.write("ERROR: - in mutation choice\n")
				continue
			log.write("CHUNK MUTATION: attempting" + str(mutations[i][0]) + "to" + mutations[i][1] + "\n")
			(position, mutation_type) = calculate_mutation_for_pose(sequence, mutations[i][0], mutations[i][1], 0)
			try:
				mutate_residue(test_pose, position, mutations[i][1])
			except:
				log.write("ERROR: Could not mutate structure \n")
				log.close()
				raise
			mutated_sequence = update_seq_string(sequence, mutations[i][1], mutations[i][0])
			log_positions.append(mutations[i][0])
			log_mutations.append(mutations[i][1])
			if (sequence[mutations[i][0]] == test_pose.sequence()[position-1]): #Pose position is 1 indexed, but strings are 0 indexed
				log.write("CHUNK MUTATION: failed mutation " + sequence + " vs. " + test_pose.sequence() + "\n")	
				return (-1, -1, -1)
			sequence = mutated_sequence
			log.write("CHUNK MUTATION: sequence =" + sequence + "\n")	
			i += 1
	else :
		#Fail state
		return (-1, -1, -1)

	#Update sequence and pose for this iteration
	mutated_sequence = sequence
	pose.assign(test_pose)
	return (mutated_sequence, log_positions, log_mutations)

def d_make_random_mutation(pose, sequence):
	"""
	Makes a random mutation
	Attempts random mutation until the structure successfully accepts one
	"""
	(position,pose_position, mutation) = random_mutation(sequence)
	#(pose_position, mutation_type) = calculate_mutation_for_pose(sequence, position, mutation, 0)
	try:
		mutate_residue(pose, pose_position, mutation)
	except:
		log.write("ERROR: Could not mutate structure \n")
		log.close()
		raise
	while(sequence[position] == (pose.sequence())[pose_position-1]):
		(position, pose_position, mutation) = random_mutation(sequence)
		#(pose_position, mutation_type) = calculate_mutation_for_pose(sequence, position, mutation, 0)
		try:
			mutate_residue(pose, pose_position, mutation)
		except:
			log.write("ERROR: Could not mutate structure \n")
			log.close()
			raise
	mutated_sequence = update_seq_string(sequence, mutation, position)
	log.write("MUTATION: " + str(position) + "," + mutation + "\n")
	return (mutated_sequence, position, mutation)


def get_current_structure():
	"""
	Retrieve the current pose and sequence or this iteration
	"""
	if (number_of_resets > 20):
		print "Max number of resets hit\n"
		return (0, 0, 0)

	curr_idx = max(0, archive_idx - 1)
	pose = Pose()
	pose.assign(pose_archive[curr_idx])
	sequence = sequence_archive[curr_idx]
	cover = cover_archive[curr_idx]
	return (pose, sequence, cover)	

#Returns the number of iterations 
def calc_convergence(curr_iter, cap_iters):
	"""
	New convergence measure: keep a sorted list of all the coverages. 
	If at any time the spread, the top ten coverage values of this iteration 
	equals the top ten of the previous iteration set a warning flag
	and update the iteration cap to terminate after 20 iterations.
	If after 20 iterations this doesn't change, then quit. 
	"""
	global same_top_ten_count
	global convergence_flag
	cover_calc = sorted(cover_archive, reverse=True)
	i = 0
	#while i < len(prev_top_ten_cover):
	while (i < len(prev_top_ten_cover)):
		print str(cover_calc[i]) + "\n"
		print str(prev_top_ten_cover[i]) + "\n"
		if (abs(cover_calc[i] - prev_top_ten_cover[i]) > 0.0001):
			convergence_flag = False
			return (cap_iters, convergence_flag)
		i += 1

	if (convergence_flag is False):
		convergence_flag = True
		return (curr_iter + 50, convergence_flag)
	else:
		return (cap_iters, convergence_flag)

def update_archives(pose, sequence, cover): 
	"""
		Update the archives for a newly accepted structure. 
		Store all the poses, sequences and coverage for each iteration.
		Store the top ten coverage values for the previous 10 iterations prior to updating
	"""
	global sequence_archive
	global pose_archive
	global cover_archive
	global archive_idx
	global reject_count
	global prev_top_ten_cover

	prev_top_ten_cover = sorted(cover_archive, reverse=True)[0:10]
	sequence_archive[archive_idx] = sequence

	pose_archive[archive_idx] = Pose()
	pose_archive[archive_idx].assign(pose)

	cover_archive[archive_idx] = cover

	#Reset reject counter
	reject_count = 0
	archive_idx += 1
	return 

def reject_archives():

	"""
	structure saving -- create a list for all potential poses.
	Save the pose from each iteration. Jump backwards when we reset and override structures
	"""	
	global archive_idx
	global reject_count
	global number_of_resets

	reject_count += 1
	prev_top_ten_cover = sorted(cover_archive, reverse=True)[0:10]

	#If five structures have been rejected, start going back in archive
	if (reject_count > 5):
		archive_idx -= 5
		archive_idx = max(0, archive_idx)
		reject_count = 0
		number_of_resets += 1
	return

def dump_intermediate_structure(pose):
	""" Dump intermediate struct to file
	"""
	global intermediate_struct_counter
	midpointFile = "output/" + name_space + "/" + str(intermediate_struct_counter) + ".pdb"
	pose.dump_pdb(midpointFile)
	intermediate_struct_counter += 1

def initialize_struct_utils(sequence, name_base, max_iters, pose, cover):
	"""
	Initialize various global values that are used by this file
	"""
	global name_space
	global sequence_archive
	global pose_archive
	global cover_archive

	name_space = name_base
	sequence_archive = [0] * max_iters
	pose_archive = [0] * max_iters
	cover_archive = [0] *max_iters

	sequence_archive[0] = sequence
	cover_archive[0] = cover
	pose_archive[0] = Pose()
	pose_archive[0].assign(pose)

def set_consensus_sequence(sequence):
	global consensus_sequence
	consensus_sequence = sequence

def set_fLog(logFile):
	global log
	log =logFile

def set_mutation_temp(temp):
	global mutation_temp
	mutation_temp = temp

def mutate_preset_sequence(pose):
	"""
	Attempts to mutate the structure into the given sequence 
	"""
	global preset_idx
	global preset_attempt
	target_aa = "null"

	#Check if consensus sequence differs from struct sequence at this position
	if (not consensus_sequence[preset_idx] == pose.sequence()[preset_idx]):

		#Attempt to mutate each position in the structure to match the given sequence
		target_aa = consensus_sequence[preset_idx]
		mutate_residue(pose, preset_idx + 1, target_aa);

	preset_idx += 1

	if (preset_idx == len(consensus_sequence) and not preset_attempt == 2):
		preset_idx = 0
		preset_attempt = 2

	if (preset_idx == len(consensus_sequence) and preset_attempt == 2):
		return (-1, -1, -1)

	return (pose, preset_idx, target_aa)


if __name__ == "__main__":
	#print calculate_mutation_for_pose("AB---CD", 5, "-")
	global sequence_master
	rosetta.init()
	# DELETION TEST
	# pose = Pose()
	# pose_from_pdb(pose, "test.pdb")
	# zero_pose = zero_pose(pose)
	# dump_pdb(zero_pose, "testZero.pdb")
	# print pose.sequence()
	# new_pose =  delete_residue(10, pose)
	# dump_pdb(new_pose, "testDel.pdb")
	# print new_pose.sequence()

	#INSERTION TEST
	# pose = Pose()
	# pose_from_pdb(pose, "test.pdb")
	# dump_pdb(pose, "testZero.pdb")
	# print pose.sequence()
	# insert_residue(pose, 10, "E")
	# dump_pdb(pose, "testIns.pdb")
	# print pose.sequence()

	#POSITION CONVERTION TEST
	# sequence_master = "SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"
	# position = 0
	# amino_acid = "-"
	# count = 0
	# print calculate_mutation_for_pose(sequence_master, position, amino_acid, count)

	#UPDATE MASTER
	# sequence_master = "SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"
	# position = 37
	# mutation = "-"
	# update_master_sequence(position, mutation)	
	# print sequence_master


	pose = Pose()
	pose_from_pdb(pose, "structures/v1v2.pdb")
	v1v2_seq = 'VKLTPLCVTLQCTNVTNNITD-------------------------------------DMRGELKN----CSFNM-T-TE--LRD-KK-QKV-YSLF-YRLDVVQINENQGNRSNNS------------------------------------------NKEYRLI---NCNTSAI-T---QA'

	# optimize_structure(pose, 100)
	# zero_pose(pose)
	# dump_pdb(pose, "../gagZero.pdb")
	# optimize_structure(pose, 100)
	# dump_pdb(pose, "../gagZeroOpt.pdb")

	# pose_from_pdb(pose, "structures/gag.pdb")
	# optimize_structure(pose, 100)
	# dump_pdb(pose, "../gagOpt.pdb")

	print pose.sequence()
	(mutated_sequence, position, mutation) = make_random_mutation(pose, v1v2_seq)
	print pose.sequence()
	print mutated_sequence
	print position
	print mutation
	print v1v2_seq[position]
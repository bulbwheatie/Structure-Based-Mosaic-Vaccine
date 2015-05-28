from rosetta import *
from utils import *
from optimizeStructure_ins import *
from PDBparser import *

possible_mutations = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# Call zero_pose() to create native structure aligned for deleteion and insertion purposes


#Keep a current [MASTER] aligned sequence (contains gaps) 
sequence_master = ""
sequence_master_archive = [0] * 5
pose_archive = [0] * 5
cover_archive = [0] * 5
intermediate_struct_counter = 0
name_space = ""
log = None
reject_count = 1

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

def make_mutation(pop, coverage_weight, mutation_length = 2):
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

	if (random.random() > point_to_chunk_prob):
		#Make a chunk mutation first
		(tmp_sequence, positions, mutations) = make_chunk_mutation(pose, sequence, cover, pop, coverage_weight, mutation_length) 
		mutation_type = "CHUNK"
		if (tmp_sequence == -1):
			#Make point mutation
			(tmp_sequence, positions, mutations) = make_point_mutation(pose, sequence, cover, coverage_weight) 
			mutation_type = "POINT"
			if (tmp_sequence == -1):
				#Make a random mutation
				(sequence, positions, mutations) = make_random_mutation(pose, sequence)
				mutation_type = "RANDOM"
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
				#Make a random mutation
				(sequence, positions, mutations) = make_random_mutation(pose, sequence)
				mutation_type = "RANDOM"
			else:
				sequence = tmp_sequence
		else:
			sequence = tmp_sequence

	#Now sequence is our mutated sequence
	#Pose is our mutated pose

	return (pose, sequence, positions, mutations, mutation_type)

def make_point_mutation(pose, sequence, cover, coverage_weight):
	(position, mutation, cover) = choose_point_mutation(sequence, cover, weight_func = coverage_weight)
	log.write("POINT MUTATION: " + str(position) + " to " + mutation + "\n")
	if (position == -1):
		return (-1, -1, -1)
	(pose_position, mutation_type) = calculate_mutation_for_pose(sequence, position, mutation, 0)
	print mutation + " at " + str(position) + "\n"
	#Check the type and make the appropriate point mutation
	if (mutation_type == "point" ):
		print "Point\n"
		mutate_residue(pose, pose_position, mutation)
		if (sequence[position] == pose.sequence()[pose_position - 1]):
			return (-1, -1, -1)
	elif (mutation_type == "insert"):
		dump_pdb(pose, "tmpRemodelPdb.pdb")
		(new_pdb, mutated_sequence) = insert_residue("tmpRemodelPdb.pdb", position, sequence)
		pose_from_pdb(pose, new_pdb)
		log.write("INSERTION\n")
		print "Insert\n"
		return (mutated_sequence, position, pose.sequence()[pose_position + 1])
	elif (mutation_type == "delete"):
		dump_pdb(pose, "tmpRemodelPdb.pdb")
		(new_pdb, mutated_sequence) = delete_residue("tmpRemodelPdb.pdb", position, sequence)
		pose_from_pdb(pose, new_pdb)
		log.write("DELETION\n")
		print "Delete\n"
		return (mutated_sequence, position, "-")

	mutated_sequence = update_seq_string(sequence, mutation, position)

	log.write("POINT MUTATION: sequence =" + str(mutated_sequence) + "\n")
	return (mutated_sequence, position, mutation)

def make_chunk_mutation(pose, sequence, cover, pop, coverage_weight, mutation_length):
	mutations = choose_n_sub_mutation(sequence, cover, pop, mut_length = mutation_length, weight_func = coverage_weight)
	log.write("CHUNK MUTATION: " + str(mutations)  + "\n")
	log_positions = []
	log_mutations = []
	print str(mutations) + "\n"
	mutated_sequence = sequence
	if (mutations is not None and mutations[0] != -1 and mutations[0][0] != -1):
		i = 0
		while i < len(mutations):
			log.write("CHUNK MUTATION: attempting" + str(mutations[i][0]) + "to" + mutations[i][1] + "\n")
			(position, mutation_type) = calculate_mutation_for_pose(sequence, mutations[i][0], mutations[i][1], 0)
			mutate_residue(pose, position, mutations[i][1])
			mutated_sequence = update_seq_string(sequence, mutations[i][1], mutations[i][0])
			log_positions.append(mutations[i][0])
			log_mutations.append(mutations[i][1])
			if (sequence[mutations[i][0]] == pose.sequence()[position-1]): #Pose position is 1 indexed, but strings are 0 indexed
				log.write("CHUNK MUTATION: failed mutation " + sequence + " vs. " + pose.sequence() + "\n")	
				return (-1, -1, -1)
			sequence = mutated_sequence
			log.write("CHUNK MUTATION: sequence =" + sequence + "\n")	
			i += 1
	else :
		#Fail state
		return (-1, -1, -1)

	#Update sequence for this iteration
	mutated_sequence = sequence
	return (mutated_sequence, log_positions, log_mutations)

def make_random_mutation(pose, sequence):
	(position,pose_position, mutation) = random_mutation(sequence)
	#(pose_position, mutation_type) = calculate_mutation_for_pose(sequence, position, mutation, 0)
	mutate_residue(pose, pose_position, mutation)
	while(sequence[position] == (pose.sequence())[pose_position-1]):
		(position, pose_position, mutation) = random_mutation(sequence)
		#(pose_position, mutation_type) = calculate_mutation_for_pose(sequence, position, mutation, 0)
		mutate_residue(pose, pose_position, mutation)

	mutated_sequence = update_seq_string(sequence, mutation, position)
	log.write("RANDOM: " + str(position) + "," + mutation + "\n")
	return (mutated_sequence, position, mutation)

def d_pose_random_mutation(testPose, pose):
	random = 0
	while (testPose.sequence() == pose.sequence()):
		random = 1
		(position, pose_position, mutation) = random_mutation(sequence_master)
		print "Random mutation " + str(position) + ", " + mutation
		(pose_position, mutation_type) = calculate_mutation_for_pose(sequence_master, position, mutation, 0)
		print "Random mutation " + str(pose_position) + ", " + mutation
		mutate_residue(pose, pose_position + 1, mutation)
	return (position, mutation, random)

def get_current_structure():
	if (pose_archive[0] != 0):
		pose = Pose()
		pose.assign(pose_archive[0])
	else:
		pose = pose_archive[0]
	sequence = sequence_master_archive[0]
	cover = cover_archive[0]
	return (pose, sequence, cover)	

def update_archives(pose, sequence, cover): 
	"""
		Update the archives for a newly accepted structure
	"""
	global sequence_master_archive
	global pose_archive
	global cover_archive
	global reject_count
	
	sequence_master_archive[1:] = sequence_master_archive[0:4]
	sequence_master_archive[0] = sequence

	pose_archive[1:] = pose_archive[0:4]
	pose_archive[0] = Pose()
	pose_archive[0].assign(pose)

	cover_archive[1:] = cover_archive[0:4]
	cover_archive[0] = cover

	#Reset reject counter
	reject_count = 0
	return 

def reject_archives():

	""" 
		Revert back one in the archives
	"""
	global sequence_master_archive
	global pose_archive
	global cover_archive
	global reject_count

	reject_count += 1

	#If three structures have been rejected, start going back in archive
	if (reject_count > 2):
		sequence_master_archive[0:4] = sequence_master_archive[1:]
		sequence_master_archive[-1] = 0
		pose_archive[0:4] = pose_archive[1:]
		pose_archive[-1] = 0
		cover_archive[0:4] = cover_archive[1:]
		cover_archive[-1] = 0
		reject_count = 0

	return
def populate_archive(pose, sequence, cover):
	global sequence_master_archive
	global pose_archive
	global cover_archive
	i = 0
	while i < len(sequence_master_archive):
		sequence_master_archive[i] = sequence
		pose_archive[i] = Pose()
		pose_archive[i].assign(pose)
		cover_archive[i] = cover
		i += 1

def d_make_mutation(pose, position, mutation, count):
	"""# Function: make_mutation()
	Determines the position and type of mutation
	Attempts to make the mutation, mutating randomly if not
	# Input: (1) Pose to mutate, (2) List of positions to mutate, (3) List of amino acids to mutate to
	# Output: Pose of mutated structure (provided pose is also changed) or None if there is an error
	# Notes: Input pose is changed if successful"""
	global sequence_master
	random = 0

	#----Validate arguments----
	if (pose is None):
		print "**ERROAR (make_mutation): Pose is null"
		return 

	if ((position) > len(sequence_master)):
		print "**ERROAR (make_mutation): Position does not exist in provided structure"
		return 

	#----Create a temporary pose
	testPose = Pose()
	testPose.assign(pose)

	log.write("DEBUG -- " + mutation + " at " + str(position) + "\n")

	#----Check what type of mutation is this (1) Insertion, (2) Deletion, (3) Point Mutation
	(pose_position, mutation_type) = calculate_mutation_for_pose(sequence_master, position, mutation, count)

	log.write("DEBUG -- " + mutation + " at " + str(position) + "\n")
	log.write("DEBUG -- " + mutation_type + " at " + str(pose_position) + "\n")
	if (mutation_type == "point" ):
		print "Point\n"
		mutate_residue(testPose, pose_position, mutation)
		pass
	elif (mutation_type == "insert"):
		insert_residue(testPose, pose_position, mutation)
		print "Insert\n"
		pass
	elif (mutation_type == "delete"):
		delete_residue(testPose, pose_position)
		print "Delete\n"

	#----Ensure mutation is successful or make a point random 
	if (testPose.sequence() == pose.sequence()):
		print "Position " + str(pose_position) + " (" + str(position) + ")\n"
		print testPose.sequence() + " vs " + pose.sequence()
		print "Making random mutaiton\n"
		(position, mutation, random) = pose_random_mutation(testPose, pose)
		mutation_type = "random"
		log.write("DEBUG -- " + mutation_type + " at " + str(position) + " (" + mutation + ")\n") 

	#----If mutation is successful, update the master sequence
	update_master_sequence(position, mutation)
	print "Master sequence = " + sequence_master + "\n"
	log.write("DEBUG -- Master sequence = " + sequence_master + "\n") 

	#----Update the provided pose and return sequence position and mutation made
	pose.assign(testPose)
	print "Assigned pose\n"
	return (position, mutation, mutation_type);

def d_su_insert_residue(pose, position, mutation):
	"""Performs an insertion into the pose right before the position
	position = position in POSE
	will insert after position specified
	"""
	position = position - 1
	sequence = pose.sequence()
	new_seq = sequence[:position] + mutation + sequence[position:]

	#Make pose from sequence
	new_pose = Pose()
	make_pose_from_sequence(new_pose, new_seq, 'fa_standard')

	#Iterate through residues, setting the angles according to original pose
	# (skipping) the deleted position
	new_idx = old_idx = 1
	while new_idx < len(new_seq):
		if new_idx == position:
			new_pose.set_phi(new_idx, 180)
			new_pose.set_psi(new_idx, 180)
			new_pose.set_omega(new_idx, 0)
			new_pose.residue(i).chi(pose.residue(i).chi())
			new_idx += 1
		else:
			new_pose.set_phi(new_idx, pose.phi(old_idx))
			new_pose.set_psi(new_idx, pose.psi(old_idx))
			new_pose.set_omega(new_idx, pose.omega(old_idx))
			new_idx += 1
			old_idx += 1
	pose.assign(new_pose)

def d_su_delete_residue(pose, position):
	""" Deletes the pose in the residue
	position = position in POSE
	"""
	sequence = pose.sequence()

	#Delete the position in the sequence
	new_seq = sequence[:position] + sequence[position + 1:]

	#Make pose from sequence
	new_pose = Pose()
	make_pose_from_sequence(new_pose, new_seq, 'fa_standard')

	#Iterate through residues, setting the angles according to original pose
	# (skipping) the deleted position
	i = 1
	while i < len(sequence):
		if i == position:
			pass
		new_pose.set_phi(i, pose.phi(i))
		new_pose.set_psi(i, pose.psi(i))
		new_pose.set_omega(i, pose.omega(i))
		new_pose.residue(i).chi(pose.residue(i).chi())
		i += 1
	pose.assign(new_pose)

def d_zero_pose(pose):
	new_pose = Pose()
	make_pose_from_sequence(new_pose, pose.sequence(), 'fa_standard')

	#Iterate through residues, setting the angles according to original pose
	# (skipping) the deleted position
	i = 1
	while i < len(pose.sequence()):
		new_pose.set_phi(i, pose.phi(i))
		new_pose.set_psi(i, pose.psi(i))
		new_pose.set_omega(i, pose.omega(i))
		new_pose.residue(i).chi(pose.residue(i).chi())
		i += 1

	#Repack the sidechains
	scorefxn = create_score_function('standard')
	packerTask = standard_packer_task(new_pose)
	packerTask.restrict_to_repacking()
	packerTask.or_include_current(True)
	packMover = PackRotamersMover(scorefxn, packerTask)
	packMover.apply(new_pose)

	pose.assign(new_pose)

def d_update_master_sequence(position, mutation):
	global sequence_master
	sequence_master = sequence_master[:position] + mutation + sequence_master[position+1:]

def get_master_sequence():
	return sequence_master

def dump_intermediate_structure(pose):
	""" Dump intermediate struct 
	"""
	global intermediate_struct_counter
	midpointFile = "output/" + name_space + "/" + str(intermediate_struct_counter) + ".pdb"
	pose.dump_pdb(midpointFile)
	intermediate_struct_counter += 1

def initialize_struct_utils(sequence, name_base):
	global sequence_master
	global name_space
	sequence_master  = sequence
	name_space = name_base


def set_fLog(logFile):
	global log
	log =logFile

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
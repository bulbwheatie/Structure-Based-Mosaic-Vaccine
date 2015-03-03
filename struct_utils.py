from rosetta import *
from utils import *

possible_mutations = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# Call zero_pose() to create native structure aligned for deleteion and insertion purposes


#Keep a current [MASTER] aligned sequence (contains gaps) 
sequence_master = ""
sequence_master_archive = ""
intermediate_struct_counter = 0
name_space = ""
log = None

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
	log.write("DEBUG -- " + amino_acid + " at " + str(position) + "\n")

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

def make_mutation(pose, position, mutation, count):
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

def insert_residue(pose, position, mutation):
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
			new_idx += 1
		else:
			new_pose.set_phi(new_idx, pose.phi(old_idx))
			new_pose.set_psi(new_idx, pose.psi(old_idx))
			new_pose.set_omega(new_idx, pose.omega(old_idx))
			new_idx += 1
			old_idx += 1
	pose.assign(new_pose)

def delete_residue(pose, position):
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
		#new_pose.set_chi(2, i, pose.chi(2, i))
		i += 1
	pose.assign(new_pose)

def zero_pose(pose):
	new_pose = Pose()
	make_pose_from_sequence(new_pose, pose.sequence(), 'fa_standard')

	#Iterate through residues, setting the angles according to original pose
	# (skipping) the deleted position
	i = 1
	while i < len(pose.sequence()):
		new_pose.set_phi(i, pose.phi(i))
		new_pose.set_psi(i, pose.psi(i))
		new_pose.set_omega(i, pose.omega(i))
		#new_pose.set_chi(2, i, pose.chi(2, i))
		i += 1

	#Repack the sidechains
	scorefxn = create_score_function('standard')
	packerTask = standard_packer_task(new_pose)
	packerTask.restrict_to_repacking()
	packerTask.or_include_current(True)
	packMover = PackRotamersMover(scorefxn, packerTask)
	packMover.apply(new_pose)

	pose.assign(new_pose)

def pose_random_mutation(testPose, pose):
	random = 0
	while (testPose.sequence() == pose.sequence()):
		random = 1
		(position, mutation) = random_mutation(sequence_master)
		print "Random mutation " + str(position) + ", " + mutation
		(pose_position, mutation_type) = calculate_mutation_for_pose(sequence_master, position, mutation, 0)
		print "Random mutation " + str(pose_position) + ", " + mutation
		mutate_residue(pose, pose_position, mutation)
	return (position, mutation, random)

def update_master_sequence(position, mutation):
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
	global sequence_master_archive
	sequence_master = sequence_master_archive = sequence
	name_space = name_base

def accept_master_sequence():
	global sequence_master_archive
	sequence_master_archive = sequence_master

def reject_master_sequence():
	global sequence_master
	sequence_master = sequence_master_archive

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
	sequence_master = "SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"
	position = 0
	amino_acid = "-"
	count = 0
	print calculate_mutation_for_pose(sequence_master, position, amino_acid, count)

	#UPDATE MASTER
	# sequence_master = "SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"
	# position = 37
	mutation = "-"
	update_master_sequence(position, mutation)	
	print sequence_master



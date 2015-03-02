from rosetta import *
from utils import *

possible_mutations = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# Call zero_pose() to create native structure aligned for deleteion and insertion purposes


#Keep a current [MASTER] aligned sequence (contains gaps) 
sequence_master = ""
intermediate_struct_counter = 0
name_space = ""

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
		(position, mutation) = random_mutation(sequence_master)

	while curr_position < position:
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

	#----Check what type of mutation is this (1) Insertion, (2) Deletion, (3) Point Mutation
	(pose_position, mutation_type) = calculate_mutation_for_pose(sequence_master, position, mutation, count)

	if (mutation_type == "point" ):
		#TODO: INSERT RESIDUE
		print "Point\n"
		mutate_residue(testPose, pose_position, mutation)
		pass
	elif (mutation_type == "insert"):
		# TODO: DELETE RESIDUE
		insert_residue(testPose, pose_position, mutation)
		print "Insert\n"
		pass
	elif (mutation_type == "delete"):
		delete_residue(testPose, pose_position)
		print "Delete\n"

	#----Ensure mutation is successful or make a point random 
	if (testPose.sequence() == pose.sequence()):
		(position, mutation, random) = pose_random_mutation(testPose, pose)
		mutation_type = "random"

	#----If mutation is successful, update the master sequence
	update_master_sequence(position, mutation)

	#----Update the provided pose and return sequence position and mutation made
	pose.assign(testPose)
	return (position, mutation, mutation_type);

def insert_residue(pose, position, mutation):
	"""Performs an insertion into the pose right after the position
	position = position in POSE
	will insert after position specified
	"""

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
	new_seq = sequence[:position-1] + sequence[position:]

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

	pose.assign(new_pose)

def pose_random_mutation(testPose, pose):
	random = 0
	while (testPose.sequence() == pose.sequence()):
		random = 1
		(position, mutation) = random_mutation(sequence_master)
		(pose_position, mutation_type) = calculate_mutation_for_pose(sequence_master, position, mutation, 0)
		mutate_residue(pose, pose_position, mutation)
	return (position, mutation, random)

def update_master_sequence(position, mutation):
	global sequence_master
	sequence_master = sequence_master[:position-1] + mutation + sequence_master[position:]


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
	sequence_master = sequence
	name_space = name_base


if __name__ == "__main__":
	#print calculate_mutation_for_pose("AB---CD", 5, "-")
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
	pose = Pose()
	pose_from_pdb(pose, "test.pdb")
	dump_pdb(pose, "testZero.pdb")
	print pose.sequence()
	insert_residue(pose, 10, "E")
	dump_pdb(pose, "testIns.pdb")
	print pose.sequence()

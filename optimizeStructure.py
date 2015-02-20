# OptimizeStructure.py 
# --------------------
# MUST CALL rosetta.init() prior to using any of these functions
#---------------------
# This file performs PyRosetta related mutations and optimizations for a protein structure 
from rosetta import *
import random 
import imp
from utils import *


possibleMutations = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def makeMutation(pose, positions, mutations):
	"""# Function: makeMutation()
	# Input: (1) Pose to mutate, (2) List of positions to mutate, (3) List of amino acids to mutate to
	# Output: Pose of mutated structure (provided pose is also changed) or None if there is an error
	# Notes: Input pose is changed if successful"""

	#Error check the parameters 
	if (not len(positions) == len(mutations)):
		print "**ERROAR (makeMutation): provide lists of equal length"
		return None

	if (pose is None):
		print "**ERROAR (makeMutation): Pose is null"
		return None

	#Create a temporary pose
	testPose = Pose()
	testPose.assign(pose)

	i = 0
	for position in positions:
		if (position > pose.total_residue()):
			print "**ERROAR (makeMutation): Position does not exist in provided structure"
			return None

		#TODO: Check if mutations[i] is valid (i.e. correct format and a natural aminio acid)

		mutate_residue(testPose, position, mutations[i])
		i +=1
	# Update the provided pose and return the mutated pose
	pose.assign(testPose)
	return pose;


def optimizeStructure(pose, positions, iters = 60):

	"""Function: optimizeStructure()
	# Input: (1) Pose to optimize, (2) List of positions that were mutated, (3) number of optimization iterations (default 60)
	# Output: Optimized structure (also changes input structure) or None if there is an error
	# Note: Attempts to optimize all the positions at once
	# TODO: Create option for attaching to PyMOL 
	"""
	#Perform error checking on parameters
	if (pose is None):
		print "**ERROAR (optimizeStructure): Pose is null"
		return None 

	for position in positions:
		if (position > pose.total_residue()):
			print "**ERROAR (optimizeStructure): Position does not exist in provided structure"
			return None

	#Create a temporary pose
	testPose = Pose()
	testPose.assign(pose)

	scorefxn = create_score_function('standard')
	kT = 1 #TODO: Put this into a constants file?
	neighborhood = 5
	numMoves = 3
	backboneAngleMax = 7
	#---------------------------------
	# CREATE THE MOVERS FOR THIS
	#---------------------------------

	mm = MoveMap()
	mm.set_bb(True)

	smallMover = SmallMover(mm, kT, numMoves)
	smallMover.angle_max('H', backboneAngleMax)
	smallMover.angle_max('E', backboneAngleMax)
	smallMover.angle_max('L', backboneAngleMax)

	shearMover = ShearMover(mm, kT, numMoves)
	shearMover.angle_max('H', backboneAngleMax)
	shearMover.angle_max('E', backboneAngleMax)
	shearMover.angle_max('L', backboneAngleMax)

	minMover = MinMover()
	minMover.movemap(mm)
	minMover.score_function(scorefxn)

	packerTask = standard_packer_task(pose)
	packerTask.restrict_to_repacking()
	packerTask.or_include_current(True)
	packMover = PackRotamersMover(scorefxn, packerTask)

	combined_mover = SequenceMover()
	combined_mover.add_mover(smallMover)
	combined_mover.add_mover(shearMover)
	combined_mover.add_mover(minMover)
	combined_mover.add_mover(packMover)

	#Use a MonteCarlo format over 60 iterations
	mc = MonteCarlo(testPose, scorefxn, kT)

	#Attach the movers to the MC object
	trial = TrialMover(combined_mover, mc)

	# Perform refinement over N number of cycles
	refinement = RepeatMover(trial, iters)

	#Performs the actual refinement process
	refinement.apply(testPose)

	mc.recover_low(testPose)
	finalScore = scorefxn(testPose)
	print finalScore
	return testPose


def ROSAIC(pdbFile, outfile, mutationGenerator, iter):
	"""# Calling ROSAIC will perform N iterations of mutation + optimization
	# INPUT: PDB file
	#		 Outfile
	#		 Mutation function that takes in a sequence and returns a list of positinos and a list of corresponding amino acids
	#
	# OUTPUT: (dumps final PDB structure)
	#         returns pdb sequence """

	rosetta.init()
	pose = Pose()
	pose_from_pdb(pose, pdbFile)
	print "Running ROSAIC"

	for i in range(0, iter):
		(positions, mutations) = mutationGenerator(pose.sequence());

		makeMutation(pose, positions, mutations);
		optimizeStructure(pose, positions)
		#TODO: Dump intermediate structures?

	pose.dump_pdb(outfile)
	return pose.sequence()

def mutationTest(sequence):
	position = [0]
	aminoAcid = [0]

	position[0] = 18
	aminoAcid[0] = possibleMutations[5]
	return (position, aminoAcid)

def mutationSelecterRandom(sequence):
	position = [0]
	aminoAcid = [0]

	all_seqs = get_all_sequences(False)
	cover = coverage(sequence, all_seqs)
	tmpCover = 0
	print cover

	while (tmpCover <= cover):
		position[0] = int(random.random() * len(sequence))
		aminoAcid[0] = possibleMutations[int(random.random() * 20)]

		tmpSequence = sequence[:position[0]] + aminoAcid[0] + sequence[(position[0] + 1):]
		tmpCover =  coverage(tmpSequence, all_seqs)
		print tmpCover

	print str(position[0]) + ", " + str(aminoAcid[0])
	return (position, aminoAcid)

if __name__ == "__main__":
    # Coverage test, should output 0.67
    print ROSAIC('nef.pdb', 'test.pdb', mutationSelecterRandom, 1)
# OptimizeStructure.py 
# --------------------
# MUST CALL rosetta.init() prior to using any of these functions
#---------------------
# This file performs PyRosetta related mutations and optimizations for a protein structure 
from rosetta import *
import random 

# Function: makeMutation()
# Input: (1) Pose to mutate, (2) List of positions to mutate, (3) List of amino acids to mutate to
# Output: Pose of mutated structure (provided pose is also changed) or None if there is an error
# Notes: Input pose is changed if successful
def makeMutation(pose, positions, mutations):

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

	for position in positions:
		if (position > pose.total_residue()):
			print "**ERROAR (makeMutation): Position does not exist in provided structure"
			return None

		#TODO: Check if mutations[i] is valid (i.e. correct format and a natural aminio acid)

		mutate_residue(testPose, positions[i], mutations[i])

	# Update the provided pose and return the mutated pose
	pose.assign(testPose)
	return pose;


# Function: optimizeStructure()
# Input: (1) Pose to optimize, (2) List of positions that were mutated, (3) number of optimization iterations (default 60)
# Output: Optimized structure (also changes input structure) or None if there is an error
# Note: Attempts to optimize all the positions at once
# TODO: Create option for attaching to PyMOL
def optimizeStructure(pose, positions, iters = 60):

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

	scorefxn = get_fa_scorefxn()
	kT = 1 #TODO: Put this into a constants file?
	neighborhood = 5
	numMoves = 3
	backboneAngleMax = 7
	# --------------------------------
	# CREATE THE MOVERS FOR THIS
	#---------------------------------

	# mmMin = MoveMap()
	# for position in positions:
	# 	#Allow perturbations for each mutation site and its 'neighborhood'
	# 	minRange = max(position - neighborhood, 1)
	# 	maxRange = min(position + neighborhood, pose.total_residue())
	# 	mmMin.set_bb_true_range(minRange, maxRange)

	# mmSmallShear = MoveMap()
	# for position in positions:
	# 	#Create a MoveMap that only allows changes at the mutation site
	# 	mmSmallShear.set_bb_true_range(position, position)

	mm = MoveMap()
	mm.set_bb(True)

	smallMover = SmallMover(mm, kT, numMoves)
	smallMover.angle_max(backboneAngleMax)

	shearMover = ShearMover(mm, kT, numMoves)
	shearMover.angle_max(backboneAngleMax)

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
	testPose.dump_pdb('./test.pdb')

	return testPose

# Calling ROSAIC will perform N iterations of mutation + optimization
# INPUT: PDB file
#		 Outfile
#		 Mutation function that takes in a sequence and returns a list of positinos and a list of corresponding amino acids
#
# OUTPUT: (dumps final PDB structure)
#         returns pdb sequence
def ROSAIC(pdbFile, outfile, mutationGenerator, iter):

	pose = pose_from_pdb(pdbFile)

	for i in range(0, iter) 
		(positions, mutations) = mutationGenerator(pose.sequence());

		makeMutation(pose, positions, mutations);
		optimizeStructure(pose, positions)

	pose.dump_pdb(outfile)
	return pose.sequence()

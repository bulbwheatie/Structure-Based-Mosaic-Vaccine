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
def optimizeStructure(pose, positions, iter = 60):
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
	maxIters = 60

	#Use a MonteCarlo format over 60 iterations
	mc = MonteCarlo(testPose, scorefxn, kT)

	for i in range(0, maxIters):
		#Perform small backbone torsion adjustments to each mutation location
		smallBackbonePerturb(testPose, positions)
		
		#Pack side chains
		repack(testPose)

		#Gradient based minimization over all locations within 5 aminio acids of a mutation site
		localMinimzation(testPose, positions)

		#Accept or reject the modified pose according to a Metropolis acceptance criterion
		mc.boltzmann()

	pose.assign(testPose)		
	return pose

# Function: localMinimzation
# Performs N iterations of gradient based local Minimization the energy function
# Optimizations are performed on the X (default 5) amino acids on either side of each mutation site
# Input: (1) pose to minimize, (2) positions of mutations, (3) scorefxn to minimize (default is Rosetta default)
def localMinimization(pose, positions, iters = 3, scorefxn = get_fa_scorefxn()):
	if (pose is None):
		print "**ERROAR (localMinimization): Pose is null"
		return None 

	#TODO: Move to some constants area
	neighborhood = 5

	#Create MinMover - performs a gradies-based minimization to the pose
	minMover = MinMover()

	#Allow perturbations for each mutation site and its 'neighborhood'
	mm = MoveMap()
	for position in positions:
		minRange = max(position - neighborhood, 0)
		maxRange = min(position + neighborhood, pose.total_residue())
		mm.set_bb_true_range(minRange, maxRange)

	minMover.movemap(mm)
	minMover.score_function(scorefxn)

	#TODO: Apply a MonteCarlo here? Does it need one? 
	for i in range(0, iters):
		minMover.appy(pose)

	return pose

# Function: repack
# Optimizes side chains from a Rotamer Library
# Input: (1) pose to repack, (2) Optional - scorefxn being used (default will be used otherwise)
def repack(pose, scorefxn = None):
	if (pose is None):
		print "**ERROAR (repack): Pose is null"
		return None 

	scorefxn = get_fa_scorefxn()
	packer_task = standard_packer_task(pose)
	packer_task.restrict_to_repacking()
	pack_mover = PackRotamersMover(scorefxn, packer_task)
	pack_mover.apply(pose)

	return pose

# Function: smallBackbonePerturb
# Performs a small or shear (80/20 chance?) backbone angle change to each of the mutation sites
# TODO: Maybe instead this should be like 5 adjustments randomly chose from the mutation sites?
def smallBackbonePerturb(pose, positions):
	if (pose is None):
		print "**ERROAR (smallBackbonePerturb): Pose is null"
		return None 

	#TODO: Standardize these somewhere in an easy to find/edit place
	kT = 1
	numMoves = 1
	maxAngle = 5

	#Create a temporary pose
	testPose = Pose()
	testPose.assign(pose)

	moveType = random.random()

	for position in positions:
		#Create a MoveMap that only allows changes at the mutation site
		mm = MoveMap()
		mm.set_bb_true_range(position)

		if (moveType > 0.8):
			#Small Mover -updates a phi or psi torsion angle by a random small angle
			small_mover = SmallMover(mm, kT, numMoves)
			small_mover.angle_max("H", maxAngle) #Not sure which of these we need
			small_mover.angle_max("E", maxAngle)
			small_mover.angle_max("L", maxAngle)
			small_mover.apply(testPose)
		else:
			#Shear Mover -updates a phi angle by a random small amount 
			# and a psi angle by the same random angle in the opposite direction
			shear_mover = ShearMover(mm, kT, numMoves)
			shear_mover.angle_max("H", maxAngle) #Not sure which of these we need
			shear_mover.angle_max("E", maxAngle)
			shear_mover.angle_max("L", maxAngle)
			shear_mover.apply(testPose)

	pose.assign(testPose)
	return pose
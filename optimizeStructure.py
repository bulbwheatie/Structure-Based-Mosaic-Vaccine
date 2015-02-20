# OptimizeStructure.py 
# --------------------
# MUST CALL rosetta.init() prior to using any of these functions
#---------------------
# This file performs PyRosetta related mutations and optimizations for a protein structure 
from rosetta import *
import random 
import imp
from utils import *
import getopt
import sys

possibleMutations = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def makeMutation(pose, position, mutation):
	"""# Function: makeMutation()
	# Input: (1) Pose to mutate, (2) List of positions to mutate, (3) List of amino acids to mutate to
	# Output: Pose of mutated structure (provided pose is also changed) or None if there is an error
	# Notes: Input pose is changed if successful"""

	if (pose is None):
		print "**ERROAR (makeMutation): Pose is null"
		return None

	#Create a temporary pose
	testPose = Pose()
	testPose.assign(pose)

	if (position > pose.total_residue()):
		print "**ERROAR (makeMutation): Position does not exist in provided structure"
		return None

	mutate_residue(testPose, position, mutation)

	# Update the provided pose and return the mutated pose
	pose.assign(testPose)
	return pose;


def optimizeStructure(pose, iters = 60):

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
	return finalScore


def ROSAIC(pdbFile, outfile, mutationGenerator, iter, logFile):
	"""# Calling ROSAIC will perform N iterations of mutation + optimization
	# INPUT: PDB file
	#		 Outfile
	#		 Mutation function that takes in a sequence and returns a list of positinos and a list of corresponding amino acids
	#
	# OUTPUT: (dumps final PDB structure)
	#         returns pdb sequence """

	print "Running ROSAIC"
	rosetta.init()
	pose = Pose()
	pose_from_pdb(pose, pdbFile)
	scorefxn = create_score_function('standard')
	print "Initial energy: " + str(scorefxn(pose))

	log = fopen(logFile, 'w')
	log.write("Initial energy: " + str(scorefxn(pose)) + "\n")

	mc = MonteCarlo(pose, scorefxn, 1)

	for i in range(0, iter):
		(position, mutation) = mutationGenerator(pose.sequence());

		makeMutation(pose, positions, mutations)
		energy = optimizeStructure(pose)
		log.write("Mutation at " + str(position) +  " to " + mutation \
			+ "; Energy = " + str(energy) + "\n")
		print "Mutation at " + str(position) +  " to " + mutation \
			+ "; Energy = " + str(energy) + "\n"
		accepted = mc.botlzmann(pose)
		if (not accepted): 
			log.write("Mutation rejected \n")
			print "Mutation rejected \n"
			
	pose.dump_pdb(outfile)
	close(log)
	return pose.sequence()

def mutationTest(sequence):
	position = 18
	aminoAcid = possibleMutations[5]
	return (position, aminoAcid)

def mutationSelecterRandom(sequence):
	position = 0
	aminoAcid = 0

	all_seqs = get_all_sequences(False)
	cover = coverage(sequence, all_seqs)
	tmpCover = 0
	print cover

	while (tmpCover <= cover):
		position = int(random.random() * len(sequence))
		aminoAcid = possibleMutations[int(random.random() * 20)]

		tmpSequence = sequence[:position] + aminoAcid + sequence[(position + 1):]
		tmpCover =  coverage(tmpSequence, all_seqs)
		print tmpCover

	print str(position) + ", " + str(aminoAcid)
	return (position, aminoAcid)

def RunROSAIC():

	pdbFile = None
	outfile = None
	logFile = None
	iters = 1

	try:
		opts, args = getopt.getopt(sys.argv[1:], "rh", ["pdbFile=", "outFile=", "iters=", "logFile"])
	except getopt.error, msg:
		print msg
		print "for help use --help"
		sys.exit(2)

	for o, a in opts:
		if o in ("--pdbFile"):
			pdbFile = a
		elif o in ("--outFile"):
			outfile = a
		elif o in ("--logFile"):
			logFile = a
		elif o in ("--iters"):
			iters = int(a);
	ROSAIC(pdbFile, outfile, mutationSelecterRandom, iters, logFile)

if __name__ == "__main__":
    print RunROSAIC()
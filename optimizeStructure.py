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
		return 

	#Create a temporary pose
	testPose = Pose()
	testPose.assign(pose)

	#Pose object is one-indexed (not zero indexed)
	if ((position + 1) > pose.total_residue()):
		print "**ERROAR (makeMutation): Position does not exist in provided structure"
		return 

	mutate_residue(testPose, position + 1, mutation)

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
		return  

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

	#Recover the lowest energy structure from the minimization process
	mc.recover_low(testPose)
	finalScore = scorefxn(testPose)

	#Modify the input pose with the minimized structure
	pose.assign(testPose)
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
	# INITIALIZE EVERYTHING
	rosetta.init()
	all_seqs = get_all_sequences(False) 
	calc_aa_ngrams(all_seqs) #Initialize ngrams for mutation chooser
	calc_pop_epitope_freq(all_seqs) #Initialize epitope freq dictionary for coverage metric

	pose = Pose()
	pose_from_pdb(pose, pdbFile)
	scorefxn = create_score_function('standard')
	print "Initial energy: " + str(scorefxn(pose))

	log = open(logFile, 'w')
	log.write("-,-,-," \
		+ str(coverage(pose.sequence())) + "," + str(scorefxn(pose)) + ",-\n")

	mc = MonteCarlo(pose, scorefxn, 100)

	for i in range(0, iter):
		#Choose a mutation
		(position, mutation, count, cover) = mutationGenerator(pose.sequence());
		if (position == -1):
			log.write("Mutations issues\n")
			pose.dump_pdb(outfile)
			log.close()
			return pose.sequence()

		#Make a mutation
		init_pose = Pose()
		init_pose.assign(pose)
		makeMutation(pose, position, mutation)

		#TODO: Check for unsuccessful mutation
		while (pose.sequence() == init_pose.sequence()):
			print "Unsuccessful mutation\n"
			pose.assign(init_pose)
			(position, aminioAcid) = random_mutation(sequence)
			makeMutation(pose, position, mutation)

		#Optimize the structure
		energy = optimizeStructure(pose)
		#Log info from this iteration
		print "Coverage " + str(cover) + " from " + str(position) +  " to " + mutation \
			+ "; Energy = " + str(energy) + "\n"
		accepted = mc.boltzmann(pose)
		print "DEBUG:: energy = " + str(scorefxn(pose)) + "\n"
		if (not accepted): 
			log.write(str(position) +  ", " + mutation \
				+ "," + str(count) + "," + str(cover) + "," + str(energy) + ",0\n")
			print "Mutation rejected \n"
		else:
			log.write(str(position) +  ", " + mutation \
				+ "," + str(count) + "," + str(cover) + "," + str(energy) + ",1\n")

	pose.dump_pdb(outfile)
	log.close()
	return pose.sequence()

def mutationTest(sequence):
	position = 18
	aminoAcid = possibleMutations[5]
	return (position, aminoAcid)

def mutationSelecter(sequence):
	"""
	Calls the choose_mutation method in utils to find a point mutation that will improve coverage
	Takes the mosaic sequence as input and outputs the position to mutate, what to mutate it to,
	number of calls to the choose mutation function to get improved coverage and what the improved coverage is
	"""
	position = 0
	aminoAcid = 0

	all_seqs = get_all_sequences(False)
	cover = coverage(sequence)
	tmpCover = cover
	print "Initial cover = " + str(cover)

	count = 0
	while (tmpCover <= cover):
		(position, aminoAcid, tmpCover) = choose_mutation(sequence, tmpCover, all_seqs)

		# No improving mutations can be found; Perform a random mutation
		if (position < 0):
			(position, aminioAcid) = random_mutation(sequence)
			tmpSequence = sequence[:position] + aminoAcid + sequence[(position + 1):]
			tmpCover = coverage(tmpSequence, all_seqs)
			return (position, aminoAcid, -1, tmpCover) #-1 count indicates a random mutation

		#Create the new sequence
		tmpSequence = sequence[:position] + aminoAcid + sequence[(position + 1):]

		#If the mutation generator doesn't compute coverage, compute it here
		if(tmpCover < 0):
			tmpCover =  coverage(tmpSequence)

		print "Cover = " + str(cover) + "; TmpCover=" + str(tmpCover)
		count +=1 
		if (count > 200):
			return (-1, aminoAcid, 200, tmpCover)

	print str(position) + ", " + str(aminoAcid)
	return (position, aminoAcid, count, tmpCover)

def RunROSAIC():
	"""
	Wrapper function to read in arguments and run ROSAIC with the current choose_mutation function from utils
	Use --pdbFile to specify the starting structure
	    --outFile to specify the name of the dumped pdb structure
	    --iters for the number of rounds of point mutations to perform
	    --logFile for the name of the log file to output data
	"""

	pdbFile = None
	outfile = None
	logFile = None
	iters = 1

	try:
		opts, args = getopt.getopt(sys.argv[1:], "rh", ["pdbFile=", "outFile=", "iters=", "logFile="])
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
	ROSAIC(pdbFile, outfile, mutationSelecter, iters, logFile)

if __name__ == "__main__":
    print RunROSAIC()
# optimize_structure.py 
# --------------------
# MUST CALL rosetta.init() prior to using any of these functions
#---------------------
# This file performs PyRosetta related mutations and optimizations for a protein structure 

from rosetta import *
import random 
import imp
from utils import *
from struct_utils import *
import getopt
import sys
import math

possible_mutations = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
intermediate_struct_counter = 0
RMSD_cutoff = 7
pop_aligned = None
pop_unaligned = None
num_mutation_choices = 2 #TODO: Wire this into struct utils
mutation_length = 2

#Constants for Rosetta's refinement process
kT = 1 #Temperature for MC
numMoves = 3 #Number of small/shear backbone torsion moves
backboneAngleMax = 7 #Maximum backbone torsion degrees that can change per iterations of refinement

energy_temp = 5
RMSD_temp = 1

def optimize_structure(pose, iters = 60):

	"""Function: optimize_structure()
	# Input: (1) Pose to optimize, (2) List of positions that were mutated, (3) number of optimization iterations (default 60)
	# Output: Optimized structure (also changes input structure) or None if there is an error
	# Note: Attempts to optimize all the positions at once
	# TODO: Create option for attaching to PyMOL 
	"""
	#Perform error checking on parameters
	if (pose is None):
		print "**ERROAR (optimize_structure): Pose is null"
		return  

	#Create a temporary pose
	testPose = Pose()
	testPose.assign(pose)

	scorefxn = create_score_function('standard')

	#---------------------------------
	# CREATE THE MOVERS FOR THIS
	#---------------------------------

	mm = MoveMap()
	mm.set_bb(True)
	#mm.set_bb_true_range(12, 16)

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

	packerTask = standard_packer_task(testPose)
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

def baseline_optimization():
	energies = []
	RMSDs = []
	pose = Pose()
	pose_from_pdb(pose, "structures/gag.pdb")
	native_pose = Pose()
	native_pose.assign(pose)
	output = open("baseline_gag.txt", 'w+')

	for iters in xrange(70):
		energy = optimize_structure(pose)
		RMSD = CA_rmsd(pose, native_pose)		
		energies.append(energy)
		RMSDs.append(RMSD)
		output.write(str(energy) + "," + str(RMSD) + "\n")

	output.close()
	return energies


def ROSAIC(pdbFile, nameBase, mutationGenerator, max_iter, sequence, coverage_weight):
	"""# Calling ROSAIC will perform N iterations of mutation + optimization
	# INPUT: PDB file
	#		 Outfile
	#		 Mutation function that takes in a sequence and returns a list of positinos and a list of corresponding amino acids
	#
	# OUTPUT: (dumps final PDB structure)
	#         returns pdb sequence
	#         dumps a intermediate PDB file whenever a random mutation is made (decrease in coverage) 
	 """

	print "Running ROSAIC\n"
	outfile = "output/" + nameBase + "/" + nameBase + ".pdb"
	logFile = "output/" + nameBase + "/" + nameBase  + ".log"
	debugFile = "output/" + nameBase + "/" + nameBase + ".debug"
	iter = max_iter

	pose = Pose() #Pose for mutation and manipulation
	native_pose = Pose() #Keep a pose of the native structure
	pose_from_pdb(native_pose, pdbFile)
	optimize_structure(native_pose)

	pose.assign(native_pose)
	scorefxn = create_score_function('standard')
	dump_intermediate_structure(pose) #Dump zero'd pose
	cover = coverage(sequence, weight_func = coverage_weight)
	print "Initial energy: " + str(scorefxn(pose))

	#Initialize struct utils
	initialize_struct_utils(sequence, nameBase, max_iter, pose, cover)

	#Keep the struct with highest coverage
	high_cover = cover
	high_pose = Pose()
	high_pose.assign(pose)
	high_seq = sequence
	high_energy = energy = scorefxn(pose)
	high_RMSD = RMSD = 0

	#Write some initial stats 
	log = open(logFile, 'w')
	debug = open(debugFile, 'w')
	log.write("Name base = " + nameBase + "\n")
	log.write("Iters = " + str(iter) + "\n")
	log.write("Soft coverage = " + str(cover) + "\n")
	log.write("Energy = " + str(scorefxn(pose)) + "\n")
	log.write("Hard Covereage = " + str(fisher_coverage(sequence, pop_aligned)) + "\n")
	log.write("Num epitopes = " + str(num_epitopes_in_mosaic(sequence, pop_aligned))+ "\n")
	log.write("Initial seq = " + sequence + "\n")
	log.write("-------START OF DATA------\n")
	set_fLog(debug)
	flag = 0

	mc = MonteCarlo(pose, scorefxn, 100)

	i= 0
	while (i < iter):
		#Choose a mutation
		# (position, mutation, count, cover) = mutationGenerator(get_master_sequence());
		# if (position == -1):
		# 	debug.write("Mutations issues\n")
		# 	pose.dump_pdb(outfile)
		# 	debug.close()
		# 	return pose.sequence()
		debug.write("ROSAIC: Initial sequence = " + str(get_current_structure)[1] + "\n")

		#Get mutations, positions, mutation types
		(pose, sequence, positions, mutations, mutation_type) = make_mutation(pop_aligned, coverage_weight, mutation_length = mutation_length)

		debug.write("ROSAIC: Mutated sequence = " + str(sequence) + "\n")

		#check for convergence based on structure
		if (pose == -1 or sequence == -1):
			log.write("------END OF DATA------\n")
			log.write("Best coverage = " + str(high_cover) + "\n")
			log.write("Best sequence = " + str(high_seq) + "\n")
			log.write("Best energy = " + str(high_energy) + "\n")
			log.write("Best RMSD = " + str(high_RMSD) + "\n")
			log.write("Hard Covereage = " + str(fisher_coverage(high_seq, pop_aligned)) + "\n")
			log.write("Num epitopes = " + str(num_epitopes_in_mosaic(high_seq, pop_aligned))+ "\n")
			high_pose.dump_pdb(outfile)
			debug.write("Structure Rejection\n")
			log.close()
			debug.close()
			return 			

		#If no mutation was found, reject the current iteration and move on 
		if (mutation_type == "NONE"):
			reject_archives()
			log.write("%.6f, %.4f, %s, %s, %s, %s, %.4f, %d, %.6f\n"%(cover, energy, str(positions), mutation_type, str(mutations), str(sequence), RMSD, 0, hard_cover))
			print "%.6f, %.4f, %s, %s, %s, %s, %.4f, %d, %.6f\n"%(cover, energy, str(positions), mutation_type, str(mutations), str(sequence), RMSD, 0, hard_cover)
			continue

		#Make a mutation
		# iter_pose = Pose() #Keep a pose (pre-mutation and optimization) for this iteration in case we need to revert
		# iter_pose.assign(pose)
		# (position, mutation, mutation_type) = make_mutation(pose, position, mutation, count) 
		# print "Mutation " + mutation + " at " + str(position)
		cover = coverage(sequence, weight_func = coverage_weight)
		if (iter %500 == 0):
			hard_cover = fisher_coverage(high_seq, pop_aligned)
		else:
			hard_cover = -1

		#Optimize the structure
		try:
			energy = optimize_structure(pose)
		except:
			debug.write("ERROR: Cannot optimize structure\n")
			log.write("------END OF DATA------\n")
			log.write("ERROR: Cannot optimize structure\n")
			print sys.exc_info()[0]
			debug.close()
			log.close()
			pose.dump_pdb(outfile)
			raise


		RMSD = CA_rmsd(pose, native_pose)
		print "Coverage " + str(cover) + "; Energy = " \
			+ str(energy) + "; RMSD =" + str(RMSD) + "\n"

		#Accept or reject the mutated structure, if rejected, revert the structure
		# (1) Check RMSD
		
		struct_accept = is_accept_struct(RMSD, scorefxn(pose), scorefxn(native_pose))
		if (struct_accept):
			#Accept the structure
			debug.write(str(cover) + "," + str(energy) + "; RMSD =" + str(RMSD) + ",1\n")
			#dump_intermediate_structure(pose) #Dump every accepted pose
			
			#Keep track of best structure
			if (cover > high_cover):
				high_cover = cover
				high_seq = sequence
				high_pose.assign(pose)
				high_energy = energy
				high_RMSD = RMSD

			 #Update the accepted pose and sequence
			update_archives(pose, sequence, cover)
		else: 
			debug.write("Rejected" + str(cover) + "," + str(energy) + "; RMSD =" + str(RMSD) + ",0\n")
			#dump_intermediate_structure(pose) #Dump every accepted pose
			#Revert the master sequence and pose -- happens auto if new struct isnt pushed
			reject_archives()
			print "Structure rejected \n"


		debug.write("ROSAIC: End of iter sequence = " + str(sequence) + "\n")
		#Output values for this iteration 
		log.write("%.6f, %.4f, %s, %s, %s, %s, %.4f, %d, %.6f\n"%(cover, energy, str(positions), mutation_type, str(mutations), str(sequence), RMSD, struct_accept, hard_cover))
		print "%.6f, %.4f, %s, %s, %s, %s, %.4f, %d, %.6f\n"%(cover, energy, str(positions), mutation_type, str(mutations), str(sequence), RMSD, struct_accept, hard_cover)

		#Check for convergence based on coverage
		# If the structure has converged, give it
		(flag, iter) = calc_convergence(i, iter)
		if (not flag):
			iter = max_iter 

		i += 1

	#Dump the final structure and return the sequence
	log.write("------END OF DATA------\n")
	log.write("Best coverage = " + str(high_cover) + "\n")
	log.write("Best sequence = " + str(high_seq) + "\n")
	log.write("Best energy = " + str(high_energy) + "\n")
	log.write("Best RMSD = " + str(high_RMSD) + "\n")
	log.write("Hard Covereage = " + str(fisher_coverage(high_seq, pop_aligned)) + "\n")
	log.write("Num epitopes = " + str(num_epitopes_in_mosaic(high_seq, pop_aligned))+ "\n")
	high_pose.dump_pdb(outfile)
	log.close()
	debug.close()
	return pose.sequence()

def is_accept_struct(RMSD, energy, native_energy):
	accept_energy = ((energy < native_energy) or (random.random() < math.exp(-(energy - native_energy/energy_temp))))
	accept_RMSD = ((RMSD < RMSD_cutoff) or (random.random() < math.exp(-(RMSD - RMSD_cutoff)/RMSD_temp)))
	return (accept_energy and accept_RMSD)

def d_mutation_selecter(sequence):
	"""
	Calls the choose_mutation method in utils to find a point mutation that will improve coverage
	Takes the mosaic sequence and name of the intermediate structure file as input. If the
	method makes a random mutation, it will first dump the structure

	Outputs the position to mutate, what to mutate it to,
	number of calls to the choose mutation function to get improved coverage and what the improved coverage is
	"""
	position = 0
	aminoAcid = 0

	cover = coverage(sequence)
	tmpCover = cover
	print "Initial cover = " + str(cover)

	count = 0
	while (tmpCover <= cover):
		(position, aminoAcid, tmpCover) = choose_mutation(sequence, tmpCover)

		# No improving mutations can be found; Tell driver to make a random mutation
		if (position < 0):
			return (0, aminoAcid, -1, tmpCover) #-1 count indicates a random mutation

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

def run_ROSAIC():
	"""
	Wrapper function to read in arguments and run ROSAIC with the current choose_mutation function from utils
	Use --pdbFile to specify the starting structure
	    --outFile to specify the name of the dumped pdb structure
	    --iters for the number of rounds of point mutations to perform
	    --logFile for the name of the log file to output data
	"""
	global pop_aligned
	global num_mutation_choices
	global pop_unaligned
	global RMSD_cutoff
	global energy_temp
	global mutation_length

	print "DEBUG -- RUN ROSAIC"
	sys.stdout.flush()

	pdbFile = None
	outfile = None
	logFile = None
	nameBase = None
	iters = 1
	fastaFile = None 
	start_i = end_i = 0
	sequence = None
	coverage_weight = None
	template = None
	mutation_temp = 0.01;

	try:
		opts, args = getopt.getopt(sys.argv[1:], "rh", ["pdbFile=", "nameBase=", "iters=", "fastaFile=", "start_i=", \
			"end_i=", "sequence=", "num_mutation_choices=", "coverage_weight=", \
			"RMSD_cutoff=", "energy_temp=", "mutation_length=", "template=", "mutation_temp="])
	except getopt.error, msg:
		print msg
		print "for help use --help"
		sys.exit(2)

	for o, a in opts:
		if o in ("--pdbFile"):
			pdbFile = a
		elif o in ("--nameBase"):
			nameBase = a
		elif o in ("--fastaFile"):
			fastaFile = a
		elif o in ("--iters"):
			iters = int(a);
		elif o in ("--start_i"):
			start_i = int(a);
		elif o in ("--end_i"):
			end_i = int(a);
		elif o in ("--energy_temp"):
			energy_temp = int(a);
		elif o in ("--mutation_length"):
			mutation_length = int(a);
		elif o in ("--num_mutation_choices"):
			num_mutation_choices = int(a);
		elif o in ("--sequence"):
			sequence = a	
		elif o in ("--template"):
			template = a			
		elif o in ("--coverage_weight"):
			coverage_weight = a
			if (coverage_weight == "exponential"):
				coverage_weight = exponential_weight
			elif (coverage_weight == "squared"):
				coverage_weight = squared_weight

		elif o in ("--RMSD_cutoff"):
			RMSD_cutoff = int(a)	 
		elif o in ("--mutation_temp"):
			mutation_temp = float(a)	

	if (coverage_weight == None):
		print "Invalid or no coverage weight specified. Defaulting to exponential_weight\n"
		coverage_weight = exponential_weight

	# INITIALIZE EVERYTHING
	rosetta.init()

	pop_aligned = read_fasta_file(fastaFile, start_i, end_i, aligned = True)

	#pop_unaligned = read_fasta_file(fastaFile, start_i, end_i, aligned = False)
	calc_pop_epitope_freq(pop_aligned)

	calc_single_freq(pop_aligned)

	set_mutation_temp(mutation_temp)
	
	#Use nameBase as the output dir for each struct
	return ROSAIC(pdbFile, nameBase, d_mutation_selecter, iters, sequence, coverage_weight)

if __name__ == "__main__":
	print "RUNNING OPTIMIZE STRUCTURE\n"
	print run_ROSAIC()
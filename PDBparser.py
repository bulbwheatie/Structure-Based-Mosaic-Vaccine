
from rosetta import *
from utils import *
import subprocess
import sys


#INDICES
i_num = 0
i_type = 1
i_alt_loc = 2
i_aa = 3
i_chain = 4
i_res_num = 5
i_x = 6
i_y = 7
i_z = 8
i_oc = 9
i_temp = 10
i_atom_type = 11

template = ""

aa_to_tri_lookup = {"A":"ALA", "C":"CYS", "D":"ASP", "E":"GLU", "F":"PHE", "G":"GLY", "H":"HIS", "I":"ILE", "K":"LYS", "L":"LEU", "M":"MET", "N":"ASN", "P":"PRO", "Q":"GLN", "R":"ARG", "S":"SER", "T":"THR", "V":"VAL", "W":"TRP", "Y":"TYR"}


def construct_amino_acid_dict(filename, dict_file):
	"""#Write a file that contains the relative off set of each atom in each amino acid  

	#For each amino acid in a file, store the type to a dictionary
	#Write values for that amino acid to file (relative to its CA)
	"""

	pdb_file = open(filename)
	outfile = open(dict_file, 'w')
	seen_aa = {}
	current_aa = ""
	curr_lines = []
	ca_x = ca_y = ca_z = 0

	#(1) Store all the lines for a residue (save CA values when we find it)
	for line in pdb_file:
		if (not line.startswith('ATOM')):
			print line
			continue
		values = parse_PDB_line(line)
		if (values[i_type] == "CA"):
			ca_x = values[i_x]
			ca_y = values[i_y]
			ca_z = values[i_z]
		if (values[i_aa] not in seen_aa):
			print "Found new aa " + values[i_aa]
			if (not(values[i_aa] != current_aa and values[i_type] == "N")):
				print "appending atom+\n"
				curr_lines.append(values)
			else: 
				print current_aa + "\n"
				print values[i_aa] + "\n"
				print values[i_type] + "\n"
				print str(curr_lines) + "\n"
				for res_line in curr_lines:
					res_line[i_x] -= ca_x
					res_line[i_y] -= ca_y
					res_line[i_z] -= ca_z
					res_line[i_num] = 0
					pdb_line = write_PDB_line(res_line)
					outfile.write(pdb_line)
					seen_aa[res_line[i_aa]] = True

				current_aa = values[i_aa]
				curr_lines = []
				curr_lines.append(values)

	
	outfile.close()
	pdb_file.close()
	return

def read_pdb_into_memory(pdb_file):
	res_dict = []

	pose = open(pdb_file)
	curr_aa = ""
	curr_res = []

	for atom in pose:
		if (atom.startswith('ATOM')):
			values = parse_PDB_line(atom)
			if (values[i_aa] != curr_aa and values[i_type] == "N"):
				if (len(curr_res) > 0):
					res_dict.append(curr_res)
				curr_aa = values[i_aa]
				curr_res = []
			curr_res.append(values)

	pose.close()
	return res_dict

def set_template_file(template_file):
	global template
	template = template_file


def insert_residue(pdb_file, seq_pos, sequence):
	"""
	Use RosettaRemodel to perform an insertion (This will modify neighboring aa as well)
	Inidivudally mutate residues to desired aa 
	Update sequence based on resulting aa's 

	INPUT:
	Pass in some temp file pdb file that will be overwritten
	Position - in the aligned sequence
	"""

	#Write a blue print file for this particular case
	pose_pos = sequence_to_pose_position(sequence, seq_pos)
	print "Inserting at " + str(pose_pos) + "\n"
	bp_file = write_blueprint_file_insert(template, pose_pos)

	#Call Rosetta Remodel
	# /home/anthill/ramos/rosetta_2014.30.57114_bundle/main/source/bin/remodel.linuxgccrelease -s gagRenum.pdb 
	#-database /home/anthill/ramos/rosetta_2014.30.57114_bundle/main/database/ -remodel:blueprint gag_blueprint.bp 
	#-remodel:quick_and_dirty -remodel::use_blueprint_sequence true -overwrite
	command = "/home/anthill/ramos/rosetta_2014.30.57114_bundle/main/source/bin/remodel.linuxgccrelease -s {struct} " + \
	"-database /home/anthill/ramos/rosetta_2014.30.57114_bundle/main/database/ -remodel:blueprint {blueprint} " + \
	"-remodel:quick_and_dirty -remodel::use_blueprint_sequence true -overwrite"
	command = command.format(struct = pdb_file, blueprint = bp_file)
	print command
	subprocess.call(command, shell=True)

	new_pdb = pdb_file.split(".")[0].split("/")[-1] + "_0001.pdb"

	#Check what the new sequence is and reflect it in the sequence file
	#Replace the mutation site in the sequence
	new_mut = aa_at_pose_position(new_pdb, pose_pos + 1)
	sequence = update_seq_string(sequence, new_mut, seq_pos)
	print "Updated sequence " +  sequence + "\n"

	#Replace aa before insertion site
	seq_pos = pose_to_sequence_position(sequence, pose_pos -1 )
	new_mut = aa_at_pose_position(new_pdb, pose_pos)
	sequence = update_seq_string(sequence, new_mut, seq_pos)
	print "Updated sequence " +  sequence + "\n"

	#Replace aa before insertion site
	seq_pos = pose_to_sequence_position(sequence, pose_pos + 1)
	new_mut = aa_at_pose_position(new_pdb, pose_pos  +2)
	sequence = update_seq_string(sequence, new_mut, seq_pos)
	print "Updated sequence " +  sequence + "\n"

	return (new_pdb, sequence)

def delete_residue(pdb_file, seq_pos, sequence):
	#Write a blue print file for this particular case
	bp_file = write_blueprint_file_delete(template, seq_pos)

	#Determine position to delete
	pose_pos = sequence_to_pose_position(sequence, seq_pos)

	#Call Rosetta Remodel
	# /home/anthill/ramos/rosetta_2014.30.57114_bundle/main/source/bin/remodel.linuxgccrelease -s gagRenum.pdb 
	#-database /home/anthill/ramos/rosetta_2014.30.57114_bundle/main/database/ -remodel:blueprint gag_blueprint.bp 
	#-remodel:quick_and_dirty -remodel::use_blueprint_sequence true -overwrite
	command = "/home/anthill/ramos/rosetta_2014.30.57114_bundle/main/source/bin/remodel.linuxgccrelease -s {struct} " + \
	"-database /home/anthill/ramos/rosetta_2014.30.57114_bundle/main/database/ -remodel:blueprint {blueprint} " + \
	"-remodel:quick_and_dirty -remodel::use_blueprint_sequence true -overwrite"
	command = command.format(struct = pdb_file, blueprint = bp_file)
	print command
	subprocess.call(command, shell=True)

	new_pdb = pdb_file.split(".")[0].split("/")[-1] + "_0001.pdb"

	#Check what the new sequence is and reflect it in the sequence file

	#Get the mutation from the mutated sequence
	sequence = update_seq_string(sequence, "-", seq_pos)
	print "Updated sequence " +  sequence + "\n"

	#Update neighbors in case they changed
	seq_pos = pose_to_sequence_position(sequence, pose_pos -1 )
	new_mut = aa_at_pose_position(new_pdb, pose_pos -1)
	sequence = update_seq_string(sequence, new_mut, seq_pos)
	print "Updated sequence " +  sequence + "\n"

	seq_pos = pose_to_sequence_position(sequence, pose_pos +1 )
	new_mut = aa_at_pose_position(new_pdb, pose_pos +1)
	sequence = update_seq_string(sequence, new_mut, seq_pos)
	print "Updated sequence " +  sequence + "\n"

	return (new_pdb, sequence)

def pose_to_sequence_position(sequence, pose_pos):
	"""
	Takes the pose position and returns the position in the sequence 

	"""
	seq_pos = 0

	for aa in sequence:
		if (aa != "-"):
			if (pose_pos<= 0):
				break;
			pose_pos -= 1
		seq_pos += 1

	return seq_pos

def sequence_to_pose_position(sequence, seq_pos):

	i =  0 
	pose_pos = 0
	while (i < seq_pos):
		if (sequence[i] != "-"):
			pose_pos += 1
		i += 1

	return pose_pos

def aa_at_pose_position(pdb_file, position):
	"""
	TESTED
	Returns the single letter amino acid at the specified position in the pose
	"""
	pose = open(pdb_file)
	curr_aa = ""
	i = 0

	for line in pose:
		if (not line.startswith('ATOM')):
			continue
		values = parse_PDB_line(line)
		if (i == position):
			return convert_to_res_aa(values[i_aa])
		if ((values[i_aa] != curr_aa and values[i_type].replace(" ", "") == "N")):
			i += 1

	return "-"

def get_pose_positions(sequence, position):
	"""
	TESTED
	Returns the mutation position and its two neighbors in the pose
	This returns zero indexed positions for the pose!!
	"""
	positions = [0,0,0]
	i = 0
	j = 0
	while (i <= position):
		if (sequence[i] != "-"):
			positions[0] = positions[1]
			positions[1] = j
			j+=1
		i += 1

	positions[2] = j 
	return positions

def get_neighbor_positions(sequence, position):
	"""
	TESTED
	Returns the unaligned positions of the mutation position and two neighbors in the sequence
	"""
	positions = [0,0,0] #INSERT First and last should be positions corresponding to aa and middle should be a blank
	i = 0
	while (i <= position):
		if (sequence[i] != "-"):
			positions[0] = positions[1]
		positions[1] = i #This position should be a blank
		i += 1

	while (i < len(sequence)):
		if (sequence[i] != "-"):
			positions[2] = i
			break
		i += 1

	return positions

def d_insert_residue(pdb_file, position, amino_acid, out_file, dict_file):
	"""
	Insert the residue in the PDB file and write a new PDB file
	Insert after position provided


	Read entire file into array (the residue) of arrays (all the atoms in the residue)
	"""
	residues = read_amino_acid_dict(dict_file)
	#print residues["ALA"]
	offset_res_num = 0 #will be 0 or 1 (based on the PDB's numbers)
	offset_atom_num = 0
	pose = read_pdb_into_memory(pdb_file)
	out_pose = open(out_file, 'w')
	find_ca_store = []
	num_atoms = 0 #Current atom that we're examining
	offset_x = 0
	offset_y = 0
	offset_z = 0
	prev_CA = []
	insert_complete = False
	pose_position = 0 #Residue position in the sequence (based on our nubmering)
	inserted_atoms = 0  #Number of atoms in inserted residue

	#Keep track of current residue so we know residue number
	curr_aa = ""

	for res in pose:
		if (pose_position == position + 1):
			#We want to insert our new residue here

			#Get the next CA entry - use this as value for inserted
			next_res = pose[pose_position]
			CA_entry = None
			N_entry = None
			C_entry = None
			N_found = False
			for atom in next_res:
				if(atom[i_type] == 'CA'):
					CA_entry = atom 
				elif(atom[i_type] == 'N' and not N_found):
					#Get the next N entry - use this as the value for the first N
					N_entry = atom
					N_found = True
				elif(atom[i_type] == 'C'):
					#Get the terminla C entry - use this as the value for the terminal C
					C_entry = atom

			#Find the previous CA entry
			prev_res = pose[pose_position-1]
			prev_CA_entry = None
			for atom in prev_res:
				if(atom[i_type] == 'CA'):
					prev_CA_entry = atom 

			#Get our new amino acid to inser
			res_to_insert = residues[convert_to_res_tri(amino_acid)]

			#Write the first N to file
			res_to_insert[0][i_x] += N_entry[i_x] 
			res_to_insert[0][i_y] += N_entry[i_y]
			res_to_insert[0][i_z] += N_entry[i_z] 
			res_to_insert[0][i_res_num] = N_entry[i_res_num]
			res_to_insert[0][i_num] = N_entry[i_num] + 1	
			out_pose.write(write_PDB_line(res_to_insert[0]))		

			for new_atom in xrange(1, len(res_to_insert)-1):
				#Write to file with offsets
				res_to_insert[new_atom][i_x] += CA_entry[i_x] 
				res_to_insert[new_atom][i_y] += CA_entry[i_y] 
				res_to_insert[new_atom][i_z] += CA_entry[i_z] 
				res_to_insert[new_atom][i_res_num] = next_res[0][i_res_num]
				res_to_insert[new_atom][i_num] = next_res[0][i_num] + new_atom
				out_pose.write(write_PDB_line(res_to_insert[new_atom]))		

			#If it's the terminal C
			res_to_insert[-1][i_x] += C_entry[i_x]
			res_to_insert[-1][i_y] += C_entry[i_y]
			res_to_insert[-1][i_z] += C_entry[i_z]
			res_to_insert[-1][i_res_num] = C_entry[i_res_num]
			res_to_insert[-1][i_num] = C_entry[i_num] + new_atom
			out_pose.write(write_PDB_line(res_to_insert[-1]))		

			#Update the offsets for all following atoms
			offset_x = (CA_entry[i_x] - float(prev_CA_entry[i_x]))  *mult
			offset_y = (CA_entry[i_y] - float(prev_CA_entry[i_y]))  *mult
			offset_z = (CA_entry[i_z] - float(prev_CA_entry[i_z]))  *mult
			offset_res_num = 1
			offset_atom_num = len(res_to_insert)

			#Mark insertion as complete
			insert_complete = True

		#Add in offsets and write to file
		for atom in res:
			atom[i_x] += offset_x
			atom[i_y] += offset_y
			atom[i_z] += offset_z
			atom[i_res_num] += offset_res_num
			atom[i_num] += offset_atom_num
			out_pose.write(write_PDB_line(atom))

		pose_position += 1

	out_pose.close()
	return

def convert_to_res_tri(aa):
	"""
	Look up the three letter code for the single amino acid request

	"""
	return aa_to_tri_lookup[aa]

def convert_to_res_aa(tri):


	for key in aa_to_tri_lookup:
		if (aa_to_tri_lookup[key] == tri):
			return key

def read_amino_acid_dict(dict_filename):
	"""
	Read in the dictionary file and store it as a dictionary
	Each residue holds an array of arrays of values for each atom
	Use indexing globals to get values
	the length of hte first array corresponds to the number of atoms in the residue
	"""
	residues = {}
	aa_dict = open(dict_filename)
	curr_aa = ""

	for line in aa_dict:
		values = parse_PDB_line(line)
		if (not(values[i_aa] != curr_aa and values[i_type] == "N")):
			residues[values[i_aa]].append(values)
		else: 
			residues[values[i_aa]] = []
			residues[values[i_aa]].append(values)
			curr_aa = values[i_aa]

	return residues

def update_xyz(values, xyz):
	"""
		Offset atom by xyz
	"""

	values[5] += xyz[0]
	values[6] += xyz[1]
	values[7] += xyz[2]
	return values

def parse_PDB_line(line):
	"""
		Parse a PDB line to return all fields
	"""
	#Get the fields in each column and then remove white space
	sys.stdout.flush()

	values = [0] * 12
	values[i_num] = int(line[6:11].replace(" ", ""))
	values[i_type] = line[12:16]
	values[i_alt_loc] = line[16]
	values[i_aa] = line[17:20].replace(" ", "")
	values[i_chain] = line[21].replace(" ", "")
	values[i_res_num] = int(line[22:26].replace(" ", ""))
	values[i_x] = float(line[30:38].replace(" ", ""))
	values[i_y] = float(line[38:46].replace(" ", ""))
	values[i_z] = float(line[46:54].replace(" ", ""))
	values[i_oc] = line[54:60].replace(" ", "")
	values[i_temp] = line[60:66].replace(" ", "")
	values[i_atom_type] = line[76:].replace(" ", "").rstrip()

	#print values
	# values[i_x] = float(values[i_x])
	# values[i_y] = float(values[i_y])
	# values[i_z] = float(values[i_z])
	# values[i_num] = int(values[i_num])
	# values[i_res_num] = int(values[i_res_num])
	return values

def renumber_PDB_file(pdb_file, out_pdb):

	pose = open(pdb_file)
	out_pose = open(out_pdb, 'w')

	curr_aa = ""
	curr_res = 0
	for line in pose:
		if (not line.startswith('ATOM')):
			out_pose.write(line)
			continue
		values = parse_PDB_line(line)
		if ((values[i_aa] != curr_aa and values[i_type].replace(" ", "") == "N")):
			curr_res += 1

		values[i_res_num] = curr_res
		out_pose.write(write_PDB_line(values))

	return out_pdb

def write_blueprint_file_insert(template, position):
	"""
	Reads in initial blueprint file and then change the position to 'x'
	"""
	temp = open(template)
	bp_filename = template + ".insert"
	bp = open(bp_filename, 'w')
	counter = 0 #position value is 0 indexed

	#Allow redesign at mutation location plus both neighbors
	for line in temp:
		if (counter == position):
			bb_type = line.split()[-1]
			bp.write("0 x " + bb_type + "\n")
			bp.write(line)
		elif (counter == position - 1):
			bp.write(line)
		else:
			#Don't remodel at any other locations
			values = line.split()
			bp.write(values[0] + " " + values[1] + " .\n")

		counter += 1

	return bp_filename

def write_blueprint_file_delete(template, position):
	"""
	Reads in initial blueprint file and then change the position to 'x'
	"""
	temp = open(template)
	bp_filename = template + ".delete"
	bp = open(template + ".delete", 'w')
	counter = 1

	#Allow for redesign at two neighboring locations
	for line in temp:

		if (counter == position):
			counter += 1
			continue
		elif (counter == position - 1):
			bp.write(line)
		elif (counter == position + 1):
			bp.write(line)
		else:
			#Don't remodel at any other locations
			values = line.split()
			bp.write(values[0] + " " + values[1] + " .\n")	
		counter +=1

	return bp_filename

def write_template_blueprint_file(pdb_file, blueprint_file):

	pose = open(pdb_file)
	blueprint = open(blueprint_file, 'w')
	curr_aa = ""

	for line in pose:
		if (not line.startswith('ATOM')):
			continue
		line = line.replace(" 1.00", " 1.00 ")
		values = parse_PDB_line(line)
		if ((values[i_aa] != curr_aa and values[i_type].replace(" ", "") == "N")):
			#Label region as loop by default
			print values
			blueprint.write(str(values[i_res_num]) + " " + convert_to_res_aa(values[i_aa]) + " L\n")

	return blueprint_file

def write_PDB_line(values):
	"""
		Format and write a PDB line to file
		i.e. ATOM     53 HD23 LEU A  58       3.396   9.937   4.220  1.00  3.72           H  
	"""
	line = "ATOM  {0:>5} {1:<4}{2:1}{3:>3} {4:1} {5:>3}    {6:>8}{7:>8}{8:>8}{9:>6}{10:>6} {11:>11}  \n"
	#pdb_line = line.format(1, 'N', 'ALA', 'A', 56, "%.3f"%5.059, "%.3f"%15.211, "%.3f"%5.394, "%.2f"%1.00, "%.2f"%s4.61, 'N')
	#print str(values)
	values[i_x] = "%.3f"%values[i_x]
	values[i_y] = "%.3f"%values[i_y]
	values[i_z] = "%.3f"%values[i_z]
	values[i_num] = str(values[i_num])
	values[i_res_num] = str(values[i_res_num])
	pdb_line = line.format(*values)
	return pdb_line


if __name__ == "__main__":
	#construct_amino_acid_dict("structures/2NEF.pdb", "data/amino_acid_dict.pdb")
	#insert_residue("structures/gag.pdb", 6, "S", "structures/insTest.pdb","data/amino_acid_dict.pdb")
	#renumber_PDB_file("structures/gag.pdb", "structures/gagRenum.pdb")
	#write_blueprint_file("structures/2ci2.renumbered.pdb", "structures/2ci2_blueprint.bp")

	#ROSETTA REMODEL COMMANDS
	# -s gagRenumber.pdb -database [path to database] -remodel:blueprint [blueprint file] \
	# -remodel:quick_and_dirty -remodel:use_pose_relax true

	# v1v2_seq = 'VKL--TPLCVT----LQCTNV-TNNITD-'
	# print len(v1v2_seq)
	# (positions) = get_neighbor_positions(v1v2_seq, 5)
	# print positions

	# positions = get_pose_positions(v1v2_seq, 5)
	# print positions

	# mutation = aa_at_pose_position("structures/gag.pdb", 3)
	# print mutation

	#RENUMBER ALL OUR PDB STRUCTURES AND GENERATE THEIR TEMPLATES
	#lines = read_pdb_into_memory("structures/gag.pdb")
	#print lines
	# renumber_PDB_file("structures/gag.pdb", "structures/gagRenum.pdb")
	# write_template_blueprint_file("structures/gagRenum.pdb", "structures/gagRenum.blueprint")

	#renumber_PDB_file("structures/2NEF.pdb", "structures/2NEFRenum.pdb")
	# write_template_blueprint_file("structures/2NEFRenum.pdb", "structures/2NEFRenum.blueprint")

	# renumber_PDB_file("structures/v1v2.pdb", "structures/v1v2Renum.pdb")
	# write_template_blueprint_file("structures/v1v2Renum.pdb", "structures/v1v2Renum.blueprint")

	# write_blueprint_file_delete("structures/gagRenum.blueprint.template", 4)
	# write_blueprint_file_insert("structures/gagRenum.blueprint.template", 30)

	#print pose_to_sequence_position("-A-BD--", 2)
	#print sequence_to_pose_position("-A--BD", 3)

	set_template_file("tmp/2NEFRenum.blueprint.template")
	pdb_file = "tmp/2NEFRenum.pdb"
	seq_pos = 2
	nef_seq = 'AWL--EA-QE-----E---E--E--VGFPVTPQVPLRPMTYKAAVDLSHFLKEKGGLEGLIHSQRRQDILDLWIYHTQGYFPDWQNYTPGPGIRYP-----------------------------------------------------------------LTFGWCYKLVPVEPEKLE-EANK---------------------------DDP-EREVLEWRFDSRLAFHHMARELHPEYF-KNA'
	# print insert_residue(pdb_file, seq_pos, nef_seq)

	print delete_residue(pdb_file, seq_pos, nef_seq)
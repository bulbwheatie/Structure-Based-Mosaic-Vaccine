import numpy as np
import matplotlib.pyplot as plt

def plot_one_run_soft_coverage_over_iters(filenames):
	for fn in filenames:
		soft_covs = []
		energies = []
		rmsds = []
		mut_pos = []
		iters, (begin_soft_cov, end_soft_cov), (begin_hard_cov, end_hard_cov), (begin_energy, end_energy), end_rmsd = read_single_log(fn, soft_coverages = soft_covs, energies = energies, rmsds = rmsds, mutation_positions = mut_pos)
		plt.plot([i for i in xrange(iters + 1)], soft_covs)
		plt.show()

def plot_avg_coverage_over_iters(filenames):
	pass

def plot_mutation_positions(filenames):
	pass

def plot_energies_over_iters(filenames):
	pass

def plot_rmsd_over_iters(filenames):
	pass

def read_single_log(filename, soft_coverages = None, energies = None, rmsds = None, mutation_positions = None):
	""" The optional arguments should be empty lists, which will be populated with values from within one run """
	raw_lines = open(filename, 'r').readlines()

	# Variables that we return
	iters = None
	begin_soft_cov = None
	end_soft_cov = None
	begin_hard_cov = None
	end_hard_cov = None
	begin_energy = None
	end_energy = None
	end_rmsd = None

	# State variables
	finished_reading_header = False
	finished_reading_data = False
	
	for l in raw_lines:
		# Header variables
		if not finished_reading_header:
			if 'Iters' in l:
				iters = int(l.split(" ")[-1])
			if 'Soft' in l:
				begin_soft_cov = float(l.split(" ")[-1])
				if soft_coverages is not None:
					soft_coverages.append(begin_soft_cov)
			if 'Hard' in l:
				begin_hard_cov = float(l.split(" ")[-1])
			if 'Energy' in l:
				begin_energy = float(l.split(" ")[-1])
				if energies is not None:
					energies.append(begin_energy)
					
		if 'START OF DATA' in l:
			finished_reading_header = True
			continue

		if finished_reading_header and not finished_reading_data:
			# Check if we've finished reading data
			if 'END OF DATA' in l:
				finished_reading_data = True
				continue
			else:
				# Read the data
				data = l.split(',')
				if soft_coverages is not None:
					soft_coverages.append(float(data[0]))
				if energies is not None:
					energies.append(float(data[1]))
				if rmsds is not None:
					rmsds.append(float(data[-2]))
				if mutation_positions is not None:
					if 'CHUNK' in l:
						temp_positions = [int(pos) for pos in data[2][2:-1].split(',')]
						mutation_positions.extend(temp_positions)
					else:
						mutation_positions.append(int(data[2]))
					
		if finished_reading_data:
			if 'Best coverage' in l:
				end_soft_cov = float(l.split(" ")[-1])
			if 'Hard' in l:
				end_hard_cov = float(l.split(" ")[-1])
			if 'energy' in l:
				end_energy = float(l.split(" ")[-1])
			if 'rmsd' in l:
				end_rmsd = float(l.split(" ")[-1])
		
	return iters, (begin_soft_cov, end_soft_cov), (begin_hard_cov, end_hard_cov), (begin_energy, end_energy), end_rmsd
			
if __name__ == "__main__":
	#print read_single_log("../anthill_output/Mar7/v1v2_6_2.2/v1v2_6_2.2.log")
	plot_one_run_soft_coverage_over_iters(["../anthill_output/Mar7/v1v2_6_2.2/v1v2_6_2.2.log"])

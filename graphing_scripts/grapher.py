import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem, entropy

def plot_one_run_soft_coverage_over_iters(filenames):
	for fn in filenames:
		soft_covs = []
		energies = []
		rmsds = []
		mut_pos = []
		iters, (begin_soft_cov, end_soft_cov), (begin_hard_cov, end_hard_cov), (begin_energy, end_energy), end_rmsd = read_single_log(fn, soft_coverages = soft_covs, energies = energies, rmsds = rmsds, mutation_positions = mut_pos)
		plt.plot([i for i in xrange(iters + 1)], soft_covs)
		plt.show()

def plot_avg_coverage_over_iters(filenames, iters, xlab, ylab, title, type = "cov"):
	avg_soft_covs = [0.0] * (iters + 1)
	avg_energies = [0.0] * (iters + 1)
	avg_rmsds = [0.0] * (iters + 1)
	all_soft_covs = []
	all_energies = []
	all_rmsds = []
	for f in filenames:
		curr_soft_covs = []
		curr_energies = []
		curr_rmsds = []
		iters, (begin_soft_cov, end_soft_cov), (begin_hard_cov, end_hard_cov), (begin_energy, end_energy), end_rmsd = read_single_log(f, soft_coverages = curr_soft_covs, energies = curr_energies, rmsds = curr_rmsds)
		for i in xrange(0, iters + 1):
			avg_soft_covs[i] += curr_soft_covs[i]
			avg_energies[i] += curr_energies[i]
			avg_rmsds[i] += curr_rmsds[i]
		all_soft_covs.append(curr_soft_covs)
		all_energies.append(curr_energies)
		all_rmsds.append(curr_rmsds)

	for i in xrange(iters + 1):
		avg_soft_covs[i] /= len(filenames)
		avg_energies[i] /= len(filenames)
		avg_rmsds[i] /= len(filenames)

	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.title(title)
	se_soft_covs = sem(np.array(all_soft_covs))
	se_energies = sem(np.array(all_energies))
	se_rmsds = sem(np.array(all_rmsds))
	if type == "cov":
		plt.errorbar([i for i in xrange(iters + 1)], avg_soft_covs, se_soft_covs, linestyle='None', marker='^')
	elif type == "energy":
		plt.errorbar([i for i in xrange(iters + 1)], avg_energies, se_energies, linestyle='None', marker='^')
	else:
		plt.errorbar([i for i in xrange(iters + 1)], avg_rmsds, se_rmsds, linestyle='None', marker='^')
	plt.show()
	return avg_soft_covs, se_soft_covs

def plot_mutation_positions(filenames, pop_file, xlab, ylab, title):
	mutation_positions = []
	for f in filenames:
		read_single_log(f, mutation_positions = mutation_positions)
	pop_seqs = open(pop_file, 'r').readlines()

	entropies = []
	for aa_i in xrange(len(pop_seqs[0])):
		counts = {}
		for seq in pop_seqs:
			if seq[aa_i] not in counts:
				counts[seq[aa_i]] = 0.0
			else:
				counts[seq[aa_i]] += 1.0 / len(pop_seqs)
		probabilities = counts.values()
		entropies.append(entropy(probabilities))

	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.title(title)
	plt.plot([i for i in xrange(len(pop_seqs[0]))], entropies, mutation_positions, [0.1] * len(mutation_positions), 'ro')
	plt.show()

def plot_energy_funcs():
	sq = [float(i ** 2) / 9 ** 2 for i in xrange(10)]
	exp = [float(2 ** i - 1) / 2 ** 9 for i in xrange(10)]
	plt.xlabel("Number of Matches")
	plt.ylabel("Score")
	plt.title("Soft Coverage Weight Functions")
	plt.plot([i for i in xrange(10)], sq, 'ro-', [i for i in xrange(10)], exp, 'b^-')
	plt.show()

def plot_cov_by_energy_func(fns_sq, soft_covs_sq, fns_exp, soft_covs_exp, xlab, ylab, title):
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.title(title)
	plt.plot()
	plt.show()


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

	if rmsds is not None:
		rmsds.append(0.0)

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
	#plot_one_run_soft_coverage_over_iters(["../anthill_output/Mar7/v1v2_6_2.2/v1v2_6_2.2.log"])
	plot_energy_funcs()
	v1v2_70iters_files_no_subs = ["../anthill_output/Mar7/v1v2_6_5.0/v1v2_6_5.0.log",
								  "../anthill_output/Mar7/v1v2_6_5.1/v1v2_6_5.1.log",
								  "../anthill_output/Mar7/v1v2_6_5.2/v1v2_6_5.2.log"]
	plot_avg_coverage_over_iters(v1v2_70iters_files_no_subs,
								 70,
								 xlab = "Iteration Count",
								 ylab = "Soft Coverage",
								 title = "v1v2 Soft Coverage vs Iteration Count (only substitutions)",
								 type = "cov")

	plot_avg_coverage_over_iters(v1v2_70iters_files_no_subs,
								 70,
								 xlab = "Iteration Count",
								 ylab = "Energy",
								 title = "v1v2 Energy vs Iteration Count (only substitutions)",
								 type = "energy")

	plot_avg_coverage_over_iters(v1v2_70iters_files_no_subs,
								 70,
								 xlab = "Iteration Count",
								 ylab = "Backbone RMSD",
								 title = "v1v2 Soft Coverage vs Iteration Count (only substitutions)",
								 type = "rmsd")

	plot_mutation_positions(v1v2_70iters_files_no_subs,
							"../logos/popv1v2.fasta",
							xlab = "Protein Position",
							ylab = "Entropy",
							title = "v1v2 Mutation Entropy and Location Choices (only substitutions)")

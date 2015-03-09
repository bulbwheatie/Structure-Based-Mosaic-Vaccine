import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem, entropy
import copy

# Same exact sequence
fisher_gag_soft = [0.8376140524388419, 0.7111610281342117, 0.8417344846369842, 0.709602563131425, 0.6699858038952113, 0.7095591644868178, 0.8417344846369842, 0.7787227258563425, 0.7694116352741763, 0.7931691330140193, 0.7156532756590396, 0.7847798937087623, 0.7704774142778548, 0.7758248459972059, 0.7327155690474859, 0.7082976932059359, 0.6195154322628946, 0.685030377384543, 0.6978039611027446, 0.6128940495100451, 0.705710481950724, 0.7348269318069599, 0.6010810490175944, 0.6548427094663221, 0.7046758901333908, 0.7153022165593454, 0.6594125509709438, 0.6003747763207559, 0.7440026882118793, 0.7719278542167767, 0.6944254928809142, 0.637302518788873, 0.6469068388418279, 0.548406899415377, 0.7371636374350352, 0.7181400228731454, 0.7168694946252179, 0.7752887613534903, 0.5701347747428602, 0.6105707342093716, 0.5696738531694705, 0.6637813261811103]

fisher_gag_hard = [0.5961862535858745, 0.3659156919341725, 0.5816520807921683, 0.38500321544788674, 0.32098337360852264, 0.3797484160141767, 0.5816520807921683, 0.4945449280279913, 0.47637219344205606, 0.5089100536247707, 0.36966408173608745, 0.5090896224041224, 0.49883376183823036, 0.4704836089881582, 0.3989840929425036, 0.4085134458760851, 0.26133471599163643, 0.30679534429062527, 0.3503672425492762, 0.26558226750166025, 0.3624093929904294, 0.4405079919288776, 0.26130435289498255, 0.31122695374457165, 0.3504931423612327, 0.39425103065102346, 0.33989160953352565, 0.24144711395192098, 0.42286584320298504, 0.4744115675694246, 0.35433593423029935, 0.29421327927040547, 0.2836717158377664, 0.1430310856269828, 0.4420416759539431, 0.39159652734223643, 0.3879135407247561, 0.5067779014643368, 0.16346576711497948, 0.25577912935133723, 0.17016211780291848, 0.33263490412377494]

# Their neffie is shorter than ours
fisher_nef_soft = [0.7054902136767102, 0.6131840246352466, 0.6144417080470943, 0.6019418970033579, 0.5715942748171305, 0.6030510216977645, 0.5851079718597971, 0.5987459843536094, 0.602082968550989, 0.47138858141209256, 0.6204677708852042, 0.567861370975658, 0.5123371262166659, 0.615880391419754, 0.6878141525232132, 0.6458334208021084, 0.5630961096239119, 0.6056535903561777, 0.5753079963400001, 0.5446235273317422, 0.6220277464087929, 0.6458334208021084, 0.5781733930910514, 0.6195576153591544, 0.5707210717911957, 0.6106652554141296, 0.5760426109049999, 0.591538967378063, 0.7066903690876097, 0.6183119588061537, 0.625462632459518, 0.6774148432762317, 0.6387617794386792, 0.6369459606884659, 0.6434679490135045, 0.5452074697123784, 0.5471387519939865, 0.5804547621645957, 0.6174426510040059, 0.6357079808538781, 0.4913412634466324, 0.5678795392092357]

fisher_nef_hard = [0.3937748329156948, 0.2994209315173829, 0.2526958320553918, 0.2601632783809904, 0.22387433225288772, 0.267394513213283, 0.18505482994801598, 0.24147163248843345, 0.23418429040143393, 0.07513001944292709, 0.27483444764899234, 0.19105398168682197, 0.14495771731912713, 0.28321288474160905, 0.37334106535837586, 0.29736667038550163, 0.19668186033036683, 0.2666026725035688, 0.22456996070548552, 0.20263918222120367, 0.2600441713941744, 0.29736667038550163, 0.20340014036672516, 0.28744917527213043, 0.21950112729113697, 0.2678150278739148, 0.2398924730307993, 0.24678753684373964, 0.4079762369922977, 0.26797570403081233, 0.2883376865060504, 0.3491889567211993, 0.3012548579996664, 0.30033489540302155, 0.3071357415496826, 0.16542099905417715, 0.17953647272992895, 0.21541243470267563, 0.2769877039500003, 0.3030284490610912, 0.11668588787053584, 0.20054660381026893]

# Their v1v2 is shorter than ours
fisher_v1v2_soft = [0.366425241006921, 0.3065077424201614, 0.3225430412283079, 0.3541868429471227, 0.31106684727174494, 0.3365903582698049, 0.35375794012879413, 0.34329224410069736, 0.32154438623388737, 0.3409512511444192, 0.25457807351097, 0.3327333932860616, 0.3332269904682607, 0.3154868018953203, 0.323886670528781, 0.27424912760524917, 0.2968987306032082, 0.2815105963209149, 0.2742080766313533, 0.3055566375708742, 0.28428352504211307, 0.2922805385157475, 0.27729640630591856, 0.33037952942218635, 0.2584389634843544, 0.2483540933993696, 0.2871331354872942, 0.3081931261758886, 0.366425241006921, 0.3515990683909957, 0.28105419308814855, 0.3429576089989608, 0.26916739594024647, 0.2788134997589447, 0.366425241006921, 0.31253281391168247, 0.34856371758098653, 0.2696226267735226, 0.33578916763521754, 0.33188409179357714, 0.28026670593473324, 0.32303616265424095]

fisher_v1v2_hard = [0.1017236195149363, 0.048859598533905924, 0.04090306473464231, 0.07761255373332476, 0.04887387594074811, 0.023852092560639526, 0.08510008154103232, 0.07075679543991449, 0.06701875268630454, 0.06429922615424041, 0.02578357624933917, 0.036821955034167284, 0.06357814892648829, 0.05605362538302116, 0.07211453171584702, 0.04169908572476565, 0.03630126283698276, 0.05986288876917835, 0.03267614580387251, 0.027909265724327736, 0.027835224390507916, 0.07308091899180852, 0.00668064467978434, 0.07276039993247956, 0.022351103513483652, 0.027250702720689003, 0.03710178209931787, 0.049512799997417756, 0.1017236195149363, 0.06855714941024776, 0.020907595973498675, 0.08788596663449652, 0.02787676178265945, 0.02816395623407098, 0.1017236195149363, 0.048000211482575175, 0.07922228322239178, 0.04056717735786433, 0.06302983776546385, 0.0694607080253324, 0.02551693005014081, 0.045639500461875565]

def plot_comparative_performance(gag_seqs, nef_seqs, v1v2_seqs):
	# Need to recalculate coverage for our runs based on ending sequence
	# and indicies that are used for nef/v1v2
	pass





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
	plt.plot([i for i in xrange(len(pop_seqs[0]))], entropies, label='entropy')
	plt.plot(mutation_positions, [0.15] * len(mutation_positions), 'ro', label='mutation')
	plt.legend()
	plt.show()

def plot_energy_funcs():
	sq = [float(i ** 2) / 9 ** 2 for i in xrange(10)]
	exp = [float(2 ** i - 1) / 2 ** 9 for i in xrange(10)]
	plt.xlabel("Number of Matches")
	plt.ylabel("Score")
	plt.title("Soft Coverage Weight Functions")
	plt.plot([i for i in xrange(10)], sq, 'ro-', label='square coverage function')
	plt.plot([i for i in xrange(10)], exp, 'b^-', label='exponential coverage function')
	plt.legend(loc='upper left')
	plt.show()

def plot_cov_by_energy_func(fns_sq, fns_exp, iters, xlab, ylab, title):
	avg_soft_covs_sq = [0.0] * (iters + 1)
	all_soft_covs_sq = []
	for f in fns_sq:
		curr_soft_covs = []
		read_single_log(f, soft_coverages = curr_soft_covs)
		for i in xrange(0, iters + 1):
			avg_soft_covs_sq[i] += curr_soft_covs[i]
		all_soft_covs_sq.append(curr_soft_covs)

	avg_soft_covs_exp = [0.0] * (iters + 1)
	all_soft_covs_exp = []
	for f in fns_exp:
		curr_soft_covs = []
		read_single_log(f, soft_coverages = curr_soft_covs)
		for i in xrange(0, iters + 1):
			avg_soft_covs_exp[i] += curr_soft_covs[i]
		all_soft_covs_exp.append(curr_soft_covs)

	for i in xrange(iters + 1):
		avg_soft_covs_sq[i] /= len(fns_sq)
		avg_soft_covs_exp[i] /= len(fns_exp)

	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.title(title)
	plt.plot([i for i in xrange(iters + 1)], avg_soft_covs_sq, 'ro', label='square coverage function')
	plt.plot([i for i in xrange(iters + 1)], avg_soft_covs_exp, 'b^', label='exponential coverage function')
	plt.legend()
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
	v1v2_70iters_files_no_subs_sq = ["../anthill_output/Mar7/v1v2_6_5.0/v1v2_6_5.0.log",
									 "../anthill_output/Mar7/v1v2_6_5.1/v1v2_6_5.1.log",
									 "../anthill_output/Mar7/v1v2_6_5.2/v1v2_6_5.2.log"]

	v1v2_70iters_files_no_subs_exp = ["../anthill_output/Mar7/v1v2_5_12.0/v1v2_5_12.0.log",
									  "../anthill_output/Mar7/v1v2_5_12.1/v1v2_5_12.1.log",
									  "../anthill_output/Mar7/v1v2_5_12.2/v1v2_5_12.2.log"]

	v1v2_files_all = copy.copy(v1v2_70iters_files_no_subs_sq)
	v1v2_files_all.extend(v1v2_70iters_files_no_subs_exp)



	# General Graphs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	plot_energy_funcs()


	# v1v2 Graphs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	plot_cov_by_energy_func(v1v2_70iters_files_no_subs_sq,
							v1v2_70iters_files_no_subs_exp,
							70,
							xlab = "Iteration",
							ylab = "Soft Coverage",
							title = "v1v2 Soft Coverage vs Energy Function")


	plot_avg_coverage_over_iters(v1v2_70iters_files_no_subs_sq,
								 70,
								 xlab = "Iteration Count",
								 ylab = "Soft Coverage",
								 title = "v1v2 Soft Coverage vs Iteration Count (only substitutions)",
								 type = "cov")

	plot_avg_coverage_over_iters(v1v2_files_all,
								 70,
								 xlab = "Iteration Count",
								 ylab = "Energy",
								 title = "v1v2 Energy vs Iteration Count (only substitutions)",
								 type = "energy")

	plot_avg_coverage_over_iters(v1v2_files_all,
								 70,
								 xlab = "Iteration Count",
								 ylab = "Backbone RMSD",
								 title = "v1v2 Soft Coverage vs Iteration Count (only substitutions)",
								 type = "rmsd")

	plot_mutation_positions(v1v2_files_all,
							"../logos/popv1v2.fasta",
							xlab = "Protein Position",
							ylab = "Entropy",
							title = "v1v2 Mutation Entropy and Location Choices (only substitutions)")
							

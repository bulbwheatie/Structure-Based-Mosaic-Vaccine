import sys, os
import utils

def score(mos, pop):
	utils.calc_pop_epitope_freq(pop)
	utils.calc_single_freq(pop)
	print "soft coverage, square"
	soft_covs = []
	for m in mos:
		soft_covs.append(utils.coverage(m, weight_func = utils.squared_weight))
	print soft_covs
	print
	print "hard coverage"
	hard_covs = []
	for m in mos:
		hard_covs.append(utils.fisher_coverage(m, pop))
	print hard_covs
	

def read_mosaic_sequences(begin, end, filename):
	lines = open(filename, 'r').read().split('>')
	mosaic_seqs = []
	for l in lines:
		if len(l.strip()) == 0:
			continue
		curr = ''.join(l.split("\n")[1:])
		mosaic_seqs.append(curr[begin:end].replace("-",""))
	return mosaic_seqs

def fisher_main(gag_pop, nef_pop, v1v2_pop, type):
	gag_mosaics = read_mosaic_sequences(280, 351, './fisher_mosaics/GAG.fasta')
	nef_mosaics = read_mosaic_sequences(23, 175, './fisher_mosaics/NEF.fasta')
	v1v2_mosaics = read_mosaic_sequences(126, 224, './fisher_mosaics/ENV.fasta')
	if type == "gag":
		print gag_mosaics[0]
		print gag_pop[0]
		print

		score(gag_mosaics, gag_pop)
	elif type == "nef":
		print nef_mosaics[0]
		print nef_pop[0]
		print

		score(nef_mosaics, nef_pop)
	else:
		print v1v2_mosaics[0]
		print v1v2_pop[0]
		print
		score(v1v2_mosaics, v1v2_pop)

if __name__ == "__main__":
	gag_pop = utils.read_fasta_file('./data/HIV-1_gag.fasta', 343, 414, aligned = False)
	nef_pop = utils.read_fasta_file('./data/HIV-1_nef.fasta', 111 + 47, 357 - 72, aligned = False)
	v1v2_pop = utils.read_fasta_file('./data/HIV-1_env.fasta', 171, 354 - 8, aligned = False)

	if sys.argv[1] == "fisher":
		fisher_main(gag_pop, nef_pop, v1v2_pop, sys.argv[2])
	else:
		# Score our own mosaics based on coverage on the restricted range
		# log file
		loc = './anthill_output/Final2/'
		temp_logfiles = os.listdir(loc)
		logfiles70 = [loc + l + '/' + l.split('/')[-1] + '.log' for l in temp_logfiles if '70' in l]
		logfiles200 = [loc + l + '/' + l.split('/')[-1] + '.log' for l in temp_logfiles if '200' in l]
		
		gag_logfiles70 = [loc + l + '/' + l.split('/')[-1] + '.log' for l in temp_logfiles if ('70' in l and 'Gag' in l)]
		nef_logfiles70 = [loc + l + '/' + l.split('/')[-1] + '.log' for l in temp_logfiles if ('70' in l and 'Nef' in l)]
		v1v2_logfiles70 = [loc + l + '/' + l.split('/')[-1] + '.log' for l in temp_logfiles if ('70' in l and 'V1V2' in l)]
		
		gag_logfiles200 = [loc + l + '/' + l.split('/')[-1] + '.log' for l in temp_logfiles if ('200' in l and 'Gag' in l)]
		nef_logfiles200 = [loc + l + '/' + l.split('/')[-1] + '.log' for l in temp_logfiles if ('200' in l and 'Nef' in l)]
		v1v2_logfiles200 = [loc + l + '/' + l.split('/')[-1] + '.log' for l in temp_logfiles if ('200' in l and 'V1V2' in l)]

		for logfile in v1v2_logfiles200:
			start_i = None
			end_i = None
			if "Gag" in logfile:
				start_i = 0
				end_i = 71
				
			elif "Nef" in logfile:
				start_i = 47
				end_i = -72

			elif "V1V2" in logfile:
				start_i = 0
				end_i = -8

			else:
				print "ERROR"
		
			logs = open(logfile, 'r').read()
			if ('Best sequence' in logs):
				our_mosaic = logs[logs.index('Best sequence'):]
				our_mosaic = our_mosaic[our_mosaic.index('=') + 2:our_mosaic.index('\n')]
				our_mosaic = our_mosaic[start_i:end_i].replace('-', '')

				print our_mosaic
				if "Gag" in logfile:
					score([our_mosaic], gag_pop)
				elif "Nef" in logfile:
					score([our_mosaic], nef_pop)
				elif "V1V2" in logfile:
					score([our_mosaic], v1v2_pop)
		

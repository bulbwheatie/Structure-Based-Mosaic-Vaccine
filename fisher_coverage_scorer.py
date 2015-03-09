import sys
import utils

def score(mos, pop):
	utils.calc_pop_epitope_freq(pop)
	utils.calc_single_freq(pop)
	print "soft coverage, exponential"
	soft_covs = []
	for m in mos:
		soft_covs.append(utils.coverage(m))
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

if __name__ == "__main__":
	gag_mosaics = read_mosaic_sequences(280, 351, './fisher_mosaics/GAG.fasta')
	nef_mosaics = read_mosaic_sequences(23, 175, './fisher_mosaics/NEF.fasta')
	v1v2_mosaics = read_mosaic_sequences(126, 224, './fisher_mosaics/ENV.fasta')
	gag_pop = utils.read_fasta_file('./data/HIV-1_gag.fasta', 343, 414, aligned = False)
	nef_pop = utils.read_fasta_file('./data/HIV-1_nef.fasta', 111 + 47, 357 - 72, aligned = False)
	v1v2_pop = utils.read_fasta_file('./data/HIV-1_env.fasta', 171, 354 - 8, aligned = False)

	if sys.argv[1] == "gag":
		score(gag_mosaics, gag_pop)
	elif sys.argv[1] == "nef":
		score(nef_mosaics, nef_pop)
	else:
		score(v1v2_mosaics, v1v2_pop)
	


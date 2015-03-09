from utils import read_fasta_file

def read_mosaic_sequences(begin, end, filename):
	lines = open(filename, 'r').read().split('>')
	mosaic_seqs = []
	for l in lines:
		print l
		raw_input()
		curr = ''.join(l.split("\n")[1:])
		mosaic_seqs.append(curr[begin:end])
	return mosaic_seqs

if __name__ == "__main__":
	gag_mosaics = read_mosaic_sequences(279, 351, './fisher_mosaics/GAG.fasta')
	nef_mosaics = read_mosaic_sequences(23, 175, './fisher_mosaics/NEF.fasta')
	v1v2_mosaics = read_mosaic_sequences(125, 224, './fisher_mosaics/ENV.fasta')
	gag_pop = read_fasta_file('./data/HIV-1_gag.fasta', 343, 414)
	nef_pop = read_fasta_file('./data/HIV-1_nef.fasta', 111, 357)
	v1v2_pop = read_fasta_file('./data/HIV-1_env.fasta', 171, 354)

	print gag_mosaics[0], gag_pop[0]
	print nef_mosaics[0], nef_pop[0]
	print v1v2_mosaics[0], v1v2_pop[0]


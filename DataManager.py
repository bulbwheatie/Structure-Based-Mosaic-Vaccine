import os

def check_incomplete(root):

	#For every dir in the specified root, check if there is an out structure
	#Print missing structures
	dirs_to_check = []
	for (tmp, dirs, files) in os.walk(root):
		dirs_to_check.append(dirs)

	found = 0
	reruns = []

	for dir in dirs_to_check[0]:
		found = 0
		for (tmp, tmp2, files) in os.walk(root + "/" + dir): 
			for fp in files:
				if (fp == dir + ".pdb"):
					found = 1
			if (not found):
				reruns.append(dir) 

	print reruns
	bashFile = open('scripts/rosaicReruns.sh', 'w')

	for run in reruns:
		(struct, iters, wf, chunk, run_num, ins) = run.split("_")
		command = "qsub scripts/{script} {name} {iters} {wf} {mut_length} \n"
		if (wf == "exp"):
			weight_func = "\"exponential\""
		else:
			weight_func = "\"squared\""
		command = command.format(script = "rosaic" + struct + "ins.sh", name = run, iters = iters, wf = weight_func, mut_length=chunk)
		bashFile.write(command)

	bashFile.close()

def write_ins_bash():

	bashFile = open("scripts/rosaicInsDelAll.sh", 'w')
	weight_func = ["\"squared\"", "\"exponential\""]
	iters= [70, 200]
	structs = ["Gag", "Nef", "V1V2"]
	run_num = [6, 7, 8, 9, 10, 11]

	for struct in structs:
		for iter in iters:
			for wf in weight_func:
				if (wf == "\"squared\""):
					wfa = "sq"
				else:
					wfa = "exp"
				for run in run_num:
					command = "qsub scripts/{script} {name} {iters} {wf} {mut_length} \n"
					command = command.format(script = "rosaic" + struct + "ins.sh", \
						name = struct + "_" + str(iter) + "_" + wfa + "_9_" + str(run) + "_i", \
						iters =iter, wf = wf, mut_length = 9)
					bashFile.write(command)
	bashFile.close()

def write_v1v2_bash():
	bashFile = open("scripts/rosaicV1V2all.sh", 'w+')
	weight_func = ["\"squared\"", "\"exponential\""]
	iters= [200]
	structs = ["V1V2"]
	run_num = [0, 1, 2, 3, 4, 5]
	chunks = [9]

	for struct in structs:
		for iter in iters:
			for wf in weight_func:
				for chunk in chunks:
					if (wf == "\"squared\""):
						wfa = "sq"
					else:
						wfa = "exp"
					for run in run_num:
						command = "qsub scripts/{script} {name} {iters} {wf} {mut_length} \n"
						command = command.format(script = "rosaic" + struct + ".sh", \
							name = struct + "_" + str(iter) + "_" + wfa + "_" + str(chunk) + "_" + str(run) + "_ni", \
							iters =iter, wf = wf, mut_length = 9)
						bashFile.write(command)

	bashFile.close()
	return
def write_v1v2_ins_bash():
	bashFile = open("scripts/rosaicV1V2insAll.sh", 'w+')
	weight_func = ["\"squared\"", "\"exponential\""]
	iters= [200]
	structs = ["V1V2"]
	run_num = [0, 1, 2, 3, 4, 5]
	chunks = [9]

	for struct in structs:
		for iter in iters:
			for wf in weight_func:
				for chunk in chunks:
					if (wf == "\"squared\""):
						wfa = "sq"
					else:
						wfa = "exp"
					for run in run_num:
						command = "qsub scripts/{script} {name} {iters} {wf} {mut_length} \n"
						command = command.format(script = "rosaic" + struct + "ins.sh", \
							name = struct + "_" + str(iter) + "_" + wfa + "_" + str(chunk) + "_" + str(run) + "_i", \
							iters =iter, wf = wf, mut_length = 9)
						bashFile.write(command)

	bashFile.close()
	return

def write_all_bash():
	bashFile = open("scripts/rosaicAll.sh", 'w+')
	weight_func = ["\"squared\"", "\"exponential\""]
	iters= [10000]
	structs = ["V1V2", "Nef"]
	run_num = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
	chunks = [6]
	RMSD_cutoffs = [2, 4]
	mutation_temp = [0.01, 0.0075, 0.005]

	for struct in structs:
		for iter in iters:
			for wf in weight_func:
				for chunk in chunks:
					for RMSD_cutoff in RMSD_cutoffs:
						if (wf == "\"squared\""):
							wfa = "sq"
						else:
							wfa = "exp"
						for run in run_num:
							for mut_temp in mutation_temp:
								command = "qsub scripts/{script} {name} {iters} {wf} {mut_length} {RMSD_cutoff} {mutation_temp}\n"
								command = command.format(script = "rosaic" + struct + ".sh", \
									name = struct + "_" + str(iter) + "_" + wfa + "_" + str(chunk) + "_" + str(run) + "_" + str(mut_temp) +"_ni", \
									iters =iter, wf = wf, mut_length = chunk, RMSD_cutoff=RMSD_cutoff, mutation_temp=mut_temp)
								bashFile.write(command)

	bashFile.close()

	return

if __name__ == "__main__":
	#check_incomplete("anthill_output/InsDel/")
	#write_ins_bash()
	write_all_bash()
	#write_v1v2_ins_bash()


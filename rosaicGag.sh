#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N rosaic
# Combining output/error messages into one file
#$ -j y
# Set memory request:
#$ -l vf=1G
# Set walltime request:
#$ -l h_rt=06:00:00
# ironfs
#$ -l ironfs
# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory /home/username
#$ -cwd
# then you tell it retain all environment variables (as the default is to scrub your environment)
#$ -V
# Now comes the command to be executed

module load python/2.7
echo loaded python...
source /home/grigoryanlab/library/PyRosetta-Release1.1-r34968.linux.64Bit/SetPyRosettaEnvironment.sh

mkdir "output/gag_5_0.1"
time python optimizeStructure.py  --pdbFile="structs/gag.pdb" --nameBase="gag_5_0.1" --iters=10 \
--fastaFile="data/HIV-1_gag.fasta" --start_i=343 --end_i=414 \
--sequence="SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"

exit 0

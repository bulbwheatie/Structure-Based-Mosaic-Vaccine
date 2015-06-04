#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N gag_rosaic
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

nameBase=$1 #i.e. gag_5_1.0
iters=$2
weight=$3
mut_length=$4
RMSD_cutoff=$5

module load python/2.7
echo loaded python...
source /home/grigoryanlab/library/PyRosetta-Release1.1-r34968.linux.64Bit/SetPyRosettaEnvironment.sh

outdir="output/gag_test"
mkdir output/gag_test
time python optimizeStructure.py  --pdbFile="structs/gag.pdb" --nameBase="gag_test" --iters=10 \
--fastaFile="data/HIV-1_gag.fasta" --start_i=343 --end_i=414 \
--sequence="SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"

exit 0

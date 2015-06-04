# Rosaic
A Structure-Driven Mosaic Vaccine Design Algorithm

Team: Ajay Kannan and Rebecca Leong

Requirements:
- Python
- PyRosetta1.1
- Rosetta3.4 or higher (paths to your Rosetta installation will need to be modified in the PDBparser.py script).

Contents:
```sh
anthill_output/ :: contains the sample output data from the V1V2 protein structure
data/ :: contains the fasta file of the aligned population sequences
fisher_mosaics/ :: mosaic sequences from the fischer's original paper for reference
graphing_scripts/ :: scripts to parse and graph data in anthill_output
logos/ :: sequence logos to show variation at each position for gag, nef and V1V2 proteins in HIV. 
scripts/ :: bash scripts to run Rosaic on the cluster (useful for sample commands)
structures/ ::pdb files of starting structures. 

DataManager.py :: generates script files to run test sets with various combinations of parameters
PDBparser.py :: (library) file to help with insertions/deletion and PDB file processing
cocktail insdel site information .txt :: identifies high importance positions to add an insertion for cocktail design
fisher_coverage_scorer.py :: python port of the coverage metric from Fischer's original paper
mutants.py :: (library) helper functions for PyRosetta substitution mutations
optimizeStructure.py :: main driver for Rosaic. Takes command line arguments and runs rosaic.
struct_utils.py :: (library) helper functions to perform mutations on a structure
utils.py :: (library) helper functions to choose mutations and calculate coverage

```

Quick start:

Download this repository, unzip the files and place them in a location accessible by the anthill environment. 

The following script will run a simple Rosaic run from an anthill environment (from the Structure-Based-Mosaic-Vaccine\ directory) that has access to the grigoryanlab PyRosetta installation.
The scripts directory includes a number of command templates that take in command line arguements. 
```sh
chmod u+x 
./scripts/rosaicSample.sh
```

Sample usage (manual/detailed use):

```sh
source [path_to_PyRosetta/SetPyRosettaEnvironment.sh]
mkdir output/gag_test
python optimizeStructure.py --pdbFile="structures/gag.pdb" --nameBase="gag_test" --iters=10 --fastaFile="data/HIV-1_env.fasta" --start_i=171 --end_i=354 --coverage_weight="squared" --mutation_length=9 --sequence="SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"
```

This command will run 10 iterations of Rosaic on the provided gag sequence and structure (these two must match). The program will output a starting structure, output structure, log file and debug file into the directory output/gag_test. The start_i and end_i correspond to the sequence's location in the aligned population. Coverage weight specifies the weighting used for fractional coverage (possible options are "exponential" [default] and "squared"). Mutation length can be any value from 1-9 and specifies the maximum possible chunk mutation size (smaller chunk mutations will also be made in the algorithm). 
Rosaic will output logs containing [soft coverage, energy, position of mutation, mutation type, mutation amino acids, sequence, RMSD, whether the structure was accepted (1 yes, 0 no), hard coverage] as a comma separated value for each iteration. Only some of the logs output hard coverage for each iteration as this signficantly slows down the program. 

Additional example runs can be taken from the bash scripts in the scripts/ directory. 

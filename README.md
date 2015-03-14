# Rosaic
A Structure-Driven Mosaic Vaccine Design Algorithm

Team: Ajay Kannan and Rebecca Leong

Requirements:
- Python
- PyRosetta1.1
- Rosetta3.4 or higher (paths to your Rosetta installation will need to be modified in the PDBparser.py script).

Rosaic is run through the optimizeStructure.py script. This driver script runs the structure algorithm as described in the paper. It will make calls to both struct_utils.py and utils.py to choose and perform mutations. The insertion and deletion variation of each script is kept in InsDel/ and is slightly less reliable. 

anthill_output/ contains data from our runs. Some runs did not run to completion due to unknown seg faults related to PyRosetta. 

data/ contains the population sequences for each structure 

InsDel/ contains the insertion-deletion variants of our main scripts. These may be less stable...

graphing_scripts/ contains the python scripts used to visualize our data log outputs.

scripts/ contains bash scripts to run jobs on a cluster

structures/ contains the input structures to Rosaic (we used gag.pdb, 2NEF.pdb and v1v2.pdb)

Sample usage (no insertion and deletions):

```sh
source [path_to_PyRosetta/SetPyRosettaEnvironment.sh]
mkdir output/gag_test
python optimizeStructure.py --pdbFile="structures/gag.pdb" --nameBase="gag_test" --iters=10 --fastaFile="data/HIV-1_env.fasta" --start_i=171 --end_i=354 --coverage_weight="squared" --mutation_length=9 --sequence="SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"
```

This command will run 10 iterations of Rosaic on the provided gag sequence and structure (these two must match). The program will output a starting structure, output structure, log file and debug file into the directory output/gag_test. The start_i and end_i correspond to the sequence's location in the aligned population. Coverage weight specifies the weighting used for fractional coverage (possible options are "exponential" [default] and "squared"). Mutation length can be any value from 1-9 and specifies the maximum possible chunk mutation size (smaller chunk mutations will also be made in the algorithm). 

Sample usage (with insertion and deletion):

```sh
source [path_to_PyRosetta]/SetPyRosettaEnvironment.sh

mkdir output/gag_test_ins
cp "gagRenum.blueprint.template" "gag_test_ins/gagRenum.blueprint.template"
cd output/gag_test_ins
python optimizeStructure_ins.py  --pdbFile="../../../structs/gagRenum.pdb" --nameBase="gag_test_ins" --iters=10 --fastaFile="../../../data/HIV-1_gag.fasta" --start_i=343 --end_i=414 --coverage_weight="squared" --mutation_length=${mut_length} --template="gagRenum.blueprint.template" --sequence="SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDSKTILKALGPGATLEEMMTACQ"
```

Insertion and deletion runs a separate script and requries blueprint files for RosettaRemodel. The program will automativally generate these from the templates, but the naming is the same, so runs should have separate directories. In order to run the inserino and deletion script, PDBparser.py must be modified to point to your installation of Rosetta. 

Rosaic will output logs containing [soft coverage, energy, position of mutation, mutation type, mutation amino acids, sequence, RMSD, whether the structure was accepted (1 yes, 0 no), hard coverage] as a comma separated value for each iteration. Only some of the logs output hard coverage for each iteration as this signficantly slows down the program. 

Additional example runs can be taken from the bash scripts in the scripts/ directory. 

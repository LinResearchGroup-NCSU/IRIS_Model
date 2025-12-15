#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Jun 10 22:48:33 2020
# File Name: cmd.preprocessing.sh
# Description: Preprocess PDB structures and generate input data for gBinder
#########################################################################

#!/bin/bash

export PDBid=$1
export protChain=$2

# Add fake CB atoms, AWSEM format requirement
cp ../../../PDBs/${PDBid}_modified.pdb native_structures_pdbs_with_virtual_cbs/native.pdb
cp ../../../PDBs/${PDBid}_Rmodified.pdb native_structures_pdbs_with_virtual_cbs/native_Rmodified.pdb
#python add_fakeCB.py  # Uncomment if fake CB atoms need to be generated

# Generate randomized sequence for the decoys
cp proteins_list.txt sequences/
cp native_structures_pdbs_with_virtual_cbs/native.pdb sequences/
cd sequences/

# Build the sequence for native.pdb
python buildseq.py native
bash cmd.cleanSequences.sh native.seq 

# RNA nomenclature modification
cp native.seq native_Rmodified.seq
bash for_gBinder_sequences.sh

# Find the indices of contacting protein-RNA residues
# Adjust cutoff for determining contacting residues (in nm)
export cutoff=1.50  # Adjusted for hydrogen bonding & stacking interactions

python find_cm_residues.py native.pdb $cutoff randomize_position_prot.txt randomize_position_RNA.txt

# Generate decoys for the RNA
rm -rf RNA_randomization
mkdir -p RNA_randomization 

cp randomize_position_RNA.txt native.seq gBinder_sequences.txt RNA_randomization/
python generate_decoy_seq_RNA.py

# Generate decoys for the protein
rm -rf prot_randomization
mkdir -p prot_randomization 

cp randomize_position_prot.txt native.seq gBinder_sequences.txt prot_randomization/
python generate_decoy_seq_prot.py

# Combine the generated RNA and protein decoys together
rm -rf CPLEX_randomization
mkdir -p CPLEX_randomization
cat RNA_randomization/native.decoys prot_randomization/native.decoys > CPLEX_randomization/native_Rmodified.decoys

cd ../

# Create the tms file, labeling RNA as '2' and protein as '1'
# Now includes expanded atom set for contact analysis
grep "CA\|CB\|NH1\|NH2\|OG\|OH\|P\|O5'\|O2'\|N1\|N3\|O2\|O6" native_structures_pdbs_with_virtual_cbs/native.pdb > tmp.txt

# Get the total number of residues
tot_resnum=$(awk 'END{print $6}' tmp.txt)

python create_tms.py sequences/RNA_randomization/randomize_position_RNA.txt $tot_resnum

# Generate and run evaluation script
gsed "s/CPLEX_NAME/$PDBid/g; s/PROT_CHAIN/$protChain/g" template_evaluate_phi.py > evaluate_phi.py
python evaluate_phi.py
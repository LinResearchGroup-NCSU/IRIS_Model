#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Jun 10 22:48:33 2020
# File Name: cmd.preprocessing.sh
# Description: 
#########################################################################
#!/bin/bash

export PDBid=$1
export protChain=$2

# Add fake CB atoms, AWSEM format requirement
cp ../../../PDBs/${PDBid}_modified.pdb native_structures_pdbs_with_virtual_cbs/native.pdb
cp ../../../PDBs/${PDBid}_Rmodified.pdb native_structures_pdbs_with_virtual_cbs/native_Rmodified.pdb
#python add_fakeCB.py

# Generate randomized sequence for the decoys;
cp proteins_list.txt sequences/
cp native_structures_pdbs_with_virtual_cbs/native.pdb sequences/
cd sequences/

# Build the sequence for native.pdb
python buildseq.py native
bash cmd.cleanSequences.sh native.seq 

# For the RNA nonmenclature modification;
cp native.seq native_Rmodified.seq

bash for_gBinder_sequences.sh

# Find the indices of contacting protein-RNA residues
# Cutoff for determining contacting residues, unit: nm
export cutoff=1.2
python find_cm_residues.py native.pdb $cutoff randomize_position_prot.txt randomize_position_RNA.txt

# Generate decoys for the RNA
rm -r RNA_randomization
mkdir -p RNA_randomization 

cp randomize_position_RNA.txt native.seq gBinder_sequences.txt RNA_randomization/

python generate_decoy_seq_RNA.py

### Generate decoys for the protein
#rm -r prot_randomization
#mkdir -p prot_randomization 
##
#cp randomize_position_prot.txt native.seq gBinder_sequences.txt prot_randomization/
##
#python generate_decoy_seq_prot.py

# Combine the generated RNA and protein decoys together
rm -r CPLEX_randomization
mkdir -p CPLEX_randomization
cat RNA_randomization/native.decoys prot_randomization/native.decoys > CPLEX_randomization/native_Rmodified.decoys

cd ../

# Create the tms file, where the RNA is labeled as '2', while the protein are labeled as '1';
grep "CA\|O5'" native_structures_pdbs_with_virtual_cbs/native.pdb > tmp.txt
# Get the total number of residues;

tot_resnum=`cat tmp.txt | awk 'END{print $6}'`
python create_tms.py sequences/RNA_randomization/randomize_position_RNA.txt $tot_resnum

gsed "s/CPLEX_NAME/$PDBid/g; s/PROT_CHAIN/$protChain/g" template_evaluate_phi.py > evaluate_phi.py
python evaluate_phi.py


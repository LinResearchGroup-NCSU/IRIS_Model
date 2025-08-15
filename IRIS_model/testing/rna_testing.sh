#!/bin/bash

export PDBid=$1

protChain="A"

# Copy the phi1_list.txt file from subdirectories
find . -mindepth 2 -name 'phi1_list.txt' -exec cp ./phi1_list.txt {} \;

# Copy proteinList.txt to relevant subdirectories
find . -type f -name "proteinList.txt" -not -path "./proteinList.txt" -exec cp ./proteinList.txt {} \;

# Copy modified PDB files
cp PDBs/${PDBid}_modified.pdb native_structures_pdbs_with_virtual_cbs/native.pdb
cp PDBs/${PDBid}_Rmodified.pdb native_structures_pdbs_with_virtual_cbs/native_Rmodified.pdb

# Uncomment this line if you need to add fake CB atoms
# python add_fakeCB.py

# Generate sequence
cp proteins_list.txt sequences/
cp native_structures_pdbs_with_virtual_cbs/native.pdb sequences/
cd sequences/

# Build the sequence for native.pdb
python buildseq.py native
bash cmd.cleanSequences.sh native.seq

# Generate modeller whole sequence
python mapDNAseq_reverse.py rna.seq rna_modeller.seq

# Combine modeller and original protein sequence
python combine_DNAPro.py

# Cutoff for determining contacting residues (in nm)
export cutoff=1.2
python find_cm_residues.py native.pdb $cutoff randomize_position_prot.txt randomize_position_RNA.txt

# Remove any previous RNA randomization and create a new folder
rm -r RNA_randomization
mkdir -p RNA_randomization

# Copy necessary files for RNA randomization
cp randomize_position_RNA.txt native.seq native.decoys RNA_randomization/

# Clean up old CPLEX randomization folder and create a new one
rm -r CPLEX_randomization
mkdir -p CPLEX_randomization

# No need to concatenate because we don't need protein-randomized sequences
cp RNA_randomization/native.decoys CPLEX_randomization/native_Rmodified.decoys

cd ../

# Create the tms file, where the RNA is labeled as '2' and the protein as '1'
grep "CA\|O5'" native_structures_pdbs_with_virtual_cbs/native_Rmodified.pdb > tmp.txt

# Get the total number of residues
tot_resnum=`cat tmp.txt | awk 'END{print $6}'`
python create_tms.py sequences/RNA_randomization/randomize_position_RNA.txt $tot_resnum

# Substitute the placeholders with actual PDBid and protChain in the template
gsed "s/CPLEX_NAME/$PDBid/g; s/PROT_CHAIN/$protChain/g" template_evaluate_phi.py > evaluate_phi.py
python evaluate_phi.py
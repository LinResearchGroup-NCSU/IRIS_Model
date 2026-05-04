###########################################################################
# This script will the residue and nucleotide IDs that are in contacts
#
# Written by Xingcheng Lin, 12/31/2022;
###########################################################################

import time
import subprocess
import os
import math
import sys
import numpy as np
import mdtraj as md
import itertools

################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
#############################################


def find_cm_residues(pdbfile, cutoff, random_position_file_protein, random_position_file_RNA):

    # Load the PDB
    pdb = md.load_pdb(pdbfile)
    # Convert into a pandas table
    pdb_table, pdb_bonds = pdb.topology.to_dataframe()

 #   print(pdb_table)

    # Select residue indices belonging to RNA
    RNA_resID = pdb_table.loc[(pdb_table['resName'] == 'A') |  (pdb_table['resName'] == 'C') | (pdb_table['resName'] == 'U') | (pdb_table['resName'] == 'G')].drop_duplicates(subset=['resSeq']).resSeq.to_numpy()
    # Select residue indices belonging to protein
    prot_resID = pdb_table.loc[(pdb_table['resName'] == 'ALA') |  (pdb_table['resName'] == 'ARG') | (pdb_table['resName'] == 'ASN') \
         | (pdb_table['resName'] == 'ASP') | (pdb_table['resName'] == 'CYS') | (pdb_table['resName'] == 'GLU') \
            | (pdb_table['resName'] == 'GLN') | (pdb_table['resName'] == 'GLY') | (pdb_table['resName'] == 'HIS') \
                | (pdb_table['resName'] == 'ILE') | (pdb_table['resName'] == 'LEU') | (pdb_table['resName'] == 'LYS') \
                    | (pdb_table['resName'] == 'MET') | (pdb_table['resName'] == 'PHE') | (pdb_table['resName'] == 'PRO') \
                        | (pdb_table['resName'] == 'SER') | (pdb_table['resName'] == 'THR') | (pdb_table['resName'] == 'TRP') \
                            | (pdb_table['resName'] == 'TYR') | (pdb_table['resName'] == 'VAL')].drop_duplicates(subset=['resSeq']).resSeq.to_numpy()
    # For mdtraj requirement, made it into 0-indexed
    RNA_resID_0_indexed = RNA_resID - 1
    prot_resID_0_indexed = prot_resID - 1

    # Select the atom indices belong to RNA
    RNA_atom_indices = pdb.topology.select("resname =~ '[AUCG]'")
    # Select the atom indices belong to protein
    prot_atom_indices = pdb.topology.select("is_protein == True")

    # Use mdtraj to calculate the contacts
    prot_RNA_respairs = list(itertools.product(RNA_resID_0_indexed, prot_resID_0_indexed))
 #   print(prot_RNA_respairs)
    prot_RNA_respair_distances, residue_pairs = md.compute_contacts(pdb, prot_RNA_respairs, scheme='closest-heavy')


    # Select those pairs that are closer than the cutoff
    prot_RNA_respair_distances = prot_RNA_respair_distances.flatten()
    cm_pair_idx = np.argwhere(prot_RNA_respair_distances < cutoff).flatten()

    np.savetxt("respair_distances.txt", prot_RNA_respair_distances, fmt="%1.3f")
    np.savetxt("res_pairs.txt", residue_pairs, fmt="%d")
    np.savetxt("cm_pair_idx.txt", cm_pair_idx, fmt="%d")

    # Recover the 1-indexed protein and RNA residue IDs that are in contact
    prot_cm_resID = []
    RNA_cm_resID = []
    for i in cm_pair_idx:
        RNA_cm_resID = np.append(RNA_cm_resID, prot_RNA_respairs[i][0] + 1)
        prot_cm_resID = np.append(prot_cm_resID, prot_RNA_respairs[i][1] + 1)

    # Remove the duplicated residue IDs
    prot_cm_resID = np.unique(prot_cm_resID).astype(int)
    prot_cm_resID = np.reshape(prot_cm_resID, (1, np.shape(prot_cm_resID)[0]))
    RNA_cm_resID = np.unique(RNA_cm_resID).astype(int)
    RNA_cm_resID = np.reshape(RNA_cm_resID, (1, np.shape(RNA_cm_resID)[0]))
#    print(prot_cm_resID)
#    print(RNA_cm_resID)

    np.savetxt(random_position_file_protein, prot_cm_resID, fmt='%d')
    np.savetxt(random_position_file_RNA, RNA_cm_resID, fmt='%d')

    return

############################################################################

if __name__ == "__main__":
    pdbfile = sys.argv[1]
    # Cutoff for determining residue contacts
    cutoff = float(sys.argv[2])
    # files for recording indicies
    random_position_file_protein = sys.argv[3]
    random_position_file_RNA = sys.argv[4]

    find_cm_residues(pdbfile, cutoff, random_position_file_protein, random_position_file_RNA)
 
    print("I slept and dreamt that life was joy. I awoke and saw that life was service. I acted and behold, service was joy.")

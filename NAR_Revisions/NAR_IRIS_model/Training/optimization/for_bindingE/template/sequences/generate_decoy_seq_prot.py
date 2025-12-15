####################################################################################
# This script generates decoy sequences based on the given native sequence.
#
# Original: Xingcheng Lin, 06/10/2020
# Cleaned + bug-fixed (gBinder skip safety, missing return protection)
####################################################################################

import os
import sys
import random
import numpy as np

sys.path.append('../../../../common_functions')
from common_function import *

###############################################################
# Utility functions
###############################################################

def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step

def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step

def RepresentsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


###############################################################
# Main decoy generation controller
###############################################################

def generate_decoy_sequences(proteins_list_file_name,
                             methods=['RNA_randomization'],
                             num_decoys=[100000],
                             randomSeed=None):

    protein_list = read_column_from_file(proteins_list_file_name, 1)

    decoys_root_directory = "./"
    os.chdir(decoys_root_directory)

    for i, method in enumerate(methods):

        if not os.path.exists(method):
            os.makedirs(method)

        os.chdir(method)

        for protein in protein_list:
            print(method, protein)

            # output file has format <protein>.decoys
            outfile_name = f"{protein}.decoys"
            output_file = open(outfile_name, 'w')

            # reproducible randomization
            random.seed(randomSeed)

            for j in range(num_decoys[i]):
                newseq = generate_decoy_sequence(protein, method=method, degree=j)

                # SAFETY FIX:
                # Sometimes newseq can be None (if gBinder conflict), so we regenerate
                # until we get a valid sequence.
                while newseq is None:
                    newseq = generate_decoy_sequence(protein, method=method, degree=j)

                output_file.write(newseq + '\n')

            output_file.close()

        os.chdir('..')


###############################################################
# Per-sequence decoy generation
###############################################################

def generate_decoy_sequence(protein, method='RNA_randomization', degree=None):

    sequences_root_directory = "../"
    seq_path = f"{sequences_root_directory}{method}/{protein}.seq"

    with open(seq_path, "r") as f:
        native_sequence = f.read().replace('\n', '')

    ###########################################################################
    # 1. Random shuffle (protein or RNA)
    ###########################################################################
    if method == 'shuffle':
        return shuffle_string(native_sequence)

    ###########################################################################
    # 2. Cyclic permutation
    ###########################################################################
    elif method == 'cyclic':
        if degree is None:
            print("Must specify degree with method cyclic")
            sys.exit()
        return cyclically_permute_string(native_sequence, degree)

    ###########################################################################
    # 3. Constrained versions
    ###########################################################################
    elif method in ('constrained_shuffle', 'constrained_cyclic'):

        native_sequence = list(native_sequence)
        tm = read_column_from_file(os.path.join(tm_root_directory, protein + '.tm'), 1)

        outside_indices = [i for i, x in enumerate(tm) if x in ('1', '3')]
        membrane_indices = [i for i, x in enumerate(tm) if x == '2']

        outside_list = get_sublist(native_sequence, outside_indices)
        membrane_list = get_sublist(native_sequence, membrane_indices)

        outside_string = ''.join(outside_list)
        membrane_string = ''.join(membrane_list)

        if method == 'constrained_shuffle':
            outside_string = shuffle_string(outside_string)
            membrane_string = shuffle_string(membrane_string)

        elif method == 'constrained_cyclic':
            if degree is None:
                print("Must specify degree with method cyclic")
                sys.exit()
            outside_string = cyclically_permute_string(outside_string, degree)
            membrane_string = cyclically_permute_string(membrane_string, degree)

        outside_list = list(outside_string)
        membrane_list = list(membrane_string)

        new_sequence = []
        for region in tm:
            if region in ('1', '3'):
                new_sequence.append(outside_list.pop(0))
            elif region == '2':
                new_sequence.append(membrane_list.pop(0))

        return ''.join(new_sequence)

    ###########################################################################
    # 4. RNA randomization (4 nucleotides)
    ###########################################################################
    elif method == 'RNA_randomization':
        sequence_toberandomized = list(native_sequence)

        resids_toberandomized = open("randomize_position_RNA.txt", 'r').readline().split()
        for resid in resids_toberandomized:
            sequence_toberandomized[int(resid) - 1] = random.choice(["a", "g", "c", "u"])

        newsequence = ''.join(sequence_toberandomized)

        # avoid good binders
        gBinder_sequences = open("gBinder_sequences.txt", 'r').read().splitlines()
        if newsequence in gBinder_sequences:
            return None
        return newsequence

    ###########################################################################
    # 5. Protein randomization (20 amino acids)
    ###########################################################################
    elif method == 'prot_randomization':
        sequence_toberandomized = list(native_sequence)

        resids_toberandomized = open("randomize_position_prot.txt", 'r').readline().split()
        aa_pool = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

        for resid in resids_toberandomized:
            sequence_toberandomized[int(resid) - 1] = random.choice(aa_pool)

        newsequence = ''.join(sequence_toberandomized)

        # avoid good binders
        gBinder_sequences = open("gBinder_sequences.txt", 'r').read().splitlines()
        if newsequence in gBinder_sequences:
            return None
        return newsequence


###############################################################
# ENTRY POINT (unchanged)
###############################################################

generate_decoy_sequences(
    "proteins_list.txt",
    methods=['prot_randomization'],
    num_decoys=[100000],
    randomSeed=0
)
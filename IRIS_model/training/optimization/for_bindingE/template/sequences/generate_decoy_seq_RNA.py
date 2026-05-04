####################################################################################
# This script will generate decoy sequences based on your given native sequence
#
# Written by Xingcheng Lin, 06/10/2020
# Revised to handle None sequences safely, 2025/12/12
####################################################################################

import math
import subprocess
import os
import time
import sys
import random
import numpy as np

sys.path.append('../../../../common_functions')
from common_function import *

################################################

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

###########################################

def generate_decoy_sequences(proteins_list_file_name, methods=['RNA_randomization'], num_decoys=[10000], randomSeed=None):
    protein_list = read_column_from_file(proteins_list_file_name, 1)
    
    decoys_root_directory = "./"
    
    os.chdir(decoys_root_directory)
    for i, method in enumerate(methods):
        if not os.path.exists(method):
            os.makedirs(method)
        os.chdir(method)

        for protein in protein_list:
            print(method, protein)
            output_file_path = "%s.decoys" % protein
            with open(output_file_path, 'w') as output_file:
                random.seed(randomSeed)
                for j in range(num_decoys[i]):
                    seq = generate_decoy_sequence(protein, method=method, degree=j)
                    if seq is not None:  # skip None sequences
                        output_file.write(seq + '\n')

        os.chdir('..')

def generate_decoy_sequence(protein, method='RNA_randomization', degree=None):

    sequences_root_directory = "../"

    seq_file_path = "%s%s.seq" % (sequences_root_directory + method + '/', protein)
    with open(seq_file_path, "r") as sequence_file:
        native_sequence = sequence_file.read().replace('\n', '')

    if method == 'shuffle':
        return shuffle_string(native_sequence)

    elif method == 'cyclic':
        if degree is None:
            print("Must specify degree with method cyclic")
            sys.exit()
        return cyclically_permute_string(native_sequence, degree)

    elif method in ['constrained_shuffle', 'constrained_cyclic']:
        native_sequence = list(native_sequence)
        tm = read_column_from_file(os.path.join(
            tm_root_directory, protein + '.tm'), 1)
        outside_indices = [i for i, x in enumerate(tm) if x == '1' or x == '3']
        outside_list = get_sublist(native_sequence, outside_indices)
        outside_string = ''.join(outside_list)
        membrane_indices = [i for i, x in enumerate(tm) if x == '2']
        membrane_list = get_sublist(native_sequence, membrane_indices)
        membrane_string = ''.join(membrane_list)
        new_sequence = []

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
        for region in tm:
            if region in ['1', '3']:
                new_sequence.append(outside_list.pop(0))
            elif region == '2':
                new_sequence.append(membrane_list.pop(0))
        return ''.join(new_sequence)

    elif method == 'RNA_randomization':
        sequence_toberandomized = list(native_sequence)
        resids_toberandomized = open("randomize_position_RNA.txt", 'r').readline().split(' ')
        for resid_toberandomized in resids_toberandomized:
            resAbbr = random.choice(list(["a", "g", "c", "u"]))
            sequence_toberandomized[int(resid_toberandomized) - 1] = resAbbr

        newsequence = ''.join(sequence_toberandomized)
        gBinder_sequences = open("gBinder_sequences.txt", 'r').read().splitlines()
        if newsequence in gBinder_sequences:
            return None  # return None if it's a gBinder
        return newsequence

    elif method == 'prot_randomization':
        sequence_toberandomized = list(native_sequence)
        resids_toberandomized = open("randomize_position_prot.txt", 'r').readline().split(' ')
        for resid_toberandomized in resids_toberandomized:
            resAbbr = random.choice(list(
                ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]))
            sequence_toberandomized[int(resid_toberandomized) - 1] = resAbbr

        newsequence = ''.join(sequence_toberandomized)
        gBinder_sequences = open("gBinder_sequences.txt", 'r').read().splitlines()
        if newsequence in gBinder_sequences:
            return None
        return newsequence

############################################

generate_decoy_sequences("proteins_list.txt", methods=['RNA_randomization'], num_decoys=[10000], randomSeed=0)
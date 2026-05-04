####################################################################################
# This script will generate decoy sequences based on your given native sequence
#
# Written by Xingcheng Lin, 06/10/2020
####################################################################################

import math
import subprocess
import os
import time
import sys

import numpy as np

sys.path.append('../../../common_functions')
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

def phi_pairwise_contact_well(res_list_tmonly, res_list_entire, neighbor_list, parameter_list, CPLEXmodeling=False, CPLEX_name='IDK'):

    r_min, r_max, kappa, min_seq_sep = parameter_list
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    phi_pairwise_contact_well = np.zeros((24, 24))
    for res1globalindex, res1 in enumerate(res_list_entire):

        res1index = get_local_index(res1)
        res1chain = get_chain(res1)

        if CPLEXmodeling:
            if (res1 in res_list_tmonly):
                for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max + 2.0):
                    res2index = get_local_index(res2)
                    res2chain = get_chain(res2)
                    res2globalindex = get_global_index(res_list_entire, res2)
                    if (CPLEX_name == '2bu1'):
                        if (res2chain == 'A'):
                            res1type = get_res_type(res_list_entire, res1)
                            res2type = get_res_type(res_list_entire, res2)
                            rij = get_interaction_distance(res1, res2)
                            phi_pairwise_contact_well[res1type][res2type] += interaction_well(
                                rij, r_min, r_max, kappa)
                            if not res1type == res2type:
                                phi_pairwise_contact_well[res2type][res1type] += interaction_well(
                                    rij, r_min, r_max, kappa)

            else:
                continue
        else:
            for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max + 2.0):
                res2index = get_local_index(res2)
                res2chain = get_chain(res2)
                res2globalindex = get_global_index(res_list_entire, res2)
                if (res1chain == res2chain and res2index - res1index >= min_seq_sep) or (res1chain != res2chain and res2globalindex > res1globalindex):
                    res1type = get_res_type(res_list_entire, res1)
                    res2type = get_res_type(res_list_entire, res2)
                    rij = get_interaction_distance(res1, res2)
                    phi_pairwise_contact_well[res1type][res2type] += interaction_well(
                        rij, r_min, r_max, kappa)
                    if not res1type == res2type:
                        phi_pairwise_contact_well[res2type][res1type] += interaction_well(
                            rij, r_min, r_max, kappa)

    phis_to_return = []
    for i in range(24):
        for j in range(i, 24):
            phis_to_return.append(phi_pairwise_contact_well[i][j])

    return phis_to_return

###########################################

def evaluate_phis_over_training_set(training_set_file, phi_list_file_name, decoy_method, tm_only=False, num_processors=1, CPLEXmodeling=False, CPLEX_name='IDK'):
    phi_list = read_phi_list(phi_list_file_name)
    print(phi_list)
    training_set = read_column_from_file(training_set_file, 1)
    print(training_set)

    # for protein in training_set:
    evaluate_phis_for_protein(training_set, phi_list, decoy_method, tm_only=tm_only, CPLEXmodeling=CPLEXmodeling, CPLEX_name=CPLEX_name)

def evaluate_phis_for_protein(training_set, phi_list, decoy_method, tm_only=False, CPLEXmodeling=False, CPLEX_name='IDK'):
    protein = training_set[0]

    print(native_structures_directory)
    structure = parse_pdb(os.path.join(native_structures_directory, protein))

    res_list_tmonly = get_res_list(structure, tm_only=True)
    res_list_entire = get_res_list(structure, tm_only=False)
    neighbor_list = get_neighbor_list(structure, tm_only=False)
    sequence = get_sequence_from_structure(structure)

    res_list_tmonly_native = res_list_tmonly
    res_list_entire_native = res_list_entire

    # Dynamically calculate max_decoys based on total decoys present in CPLEX_randomization folder
    decoy_file_path = os.path.join(decoys_root_directory, f"{decoy_method}/{protein}.decoys")
    if os.path.exists(decoy_file_path):
        decoy_sequences = read_decoy_sequences(decoy_file_path)
        max_decoys = len(decoy_sequences)
    else:
        max_decoys = 0

    for phi, parameters in phi_list:

        phi = globals()[phi]
        parameters_string = get_parameters_string(parameters)

        # Write native phi
        output_file = open(os.path.join(phis_directory, f"{phi.__name__}_{protein}_native_{parameters_string}"), 'w')
        phis_to_write = phi(res_list_tmonly_native, res_list_entire_native,
                            neighbor_list, parameters, CPLEXmodeling=CPLEXmodeling, CPLEX_name=CPLEX_name)
        output_file.write(str(phis_to_write).strip('[]').replace(',', '') + '\n')
        output_file.close()

        # Write decoy phis
        output_file = open(os.path.join(phis_directory, f"{phi.__name__}_{protein}_decoys_{decoy_method}_{parameters_string}"), 'w')
        for i_decoy, decoy_sequence in enumerate(decoy_sequences):
            if i_decoy >= max_decoys:
                break
            mutate_whole_sequence(res_list_entire, decoy_sequence)
            phis_to_write = phi(res_list_tmonly, res_list_entire,
                                neighbor_list, parameters, CPLEXmodeling=CPLEXmodeling, CPLEX_name=CPLEX_name)
            output_file.write(str(phis_to_write).strip('[]').replace(',', ' ') + '\n')
        output_file.close()

############################################

native_structures_directory = "./native_structures_pdbs_with_virtual_cbs/"
phis_directory = "./phis/"
decoys_root_directory = "./sequences/"

# Call script, max_decoys will be determined automatically based on available decoys
evaluate_phis_over_training_set(
    "proteins_list.txt", "phi1_list.txt",
    decoy_method='CPLEX_randomization', 
    tm_only=False, num_processors=1, CPLEXmodeling=True, CPLEX_name='2bu1'
)
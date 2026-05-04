####################################################################################
# Memory-safe gamma solver (auto decoy count)
# - Reads number of decoys from the phi decoy file
# - Matches decoy file to protein AND phi parameters (robust matching)
####################################################################################

import os
import sys
import numpy as np

# -------------------------------
# Correct import path (fixed)
# -------------------------------
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_dir, '../../common_functions'))
from common_function import (
    read_phi_list,
    read_native_phi,
    get_total_phis_and_parameter_string,
    read_column_from_file
)

# -------------------------------
# User settings
# -------------------------------
gammas_directory = "./gammas/randomized_decoy/"
phis_directory   = "./phis/"

manual_cutoff = 25
use_dynamic_cutoff = False
batch_size = 2000
noise_iterations = 8
relative_error_threshold = 0.06

# Ensure directories exist
for d in [gammas_directory, phis_directory]:
    os.makedirs(d, exist_ok=True)

# -------------------------------
# AUTO-DETECT DECOY FILE BASED ON PROTEIN + φ PARAMETERS (robust)
# -------------------------------
def normalize_param_string(full_name):
    """
    Given full_name like:
      'phi_pairwise_contact_well-12.0_12.0_0.7_10'
    produce a set of parameter-only variants to search for in filenames:
      ['12.0_12.0_0.7_10', '12.0-12.0-0.7-10', '-12.0_12.0_0.7_10', ...]
    """
    # remove leading "phi_pairwise_contact_well" if present
    name = full_name
    if name.startswith("phi_pairwise_contact_well"):
        name = name[len("phi_pairwise_contact_well"):]

    # strip leading separators
    name = name.lstrip("_-")

    # produce variants
    variants = set()
    variants.add(name)  # as-is, e.g. 12.0_12.0_0.7_10
    variants.add(name.replace("-", "_"))
    variants.add(name.replace("_", "-"))
    variants.add("_" + name)
    variants.add("-" + name)
    # also without decimal trailing zeros differences (in case)
    variants.add(name.replace(".0", ""))           # 12_12_0.7_10
    variants.add(name.replace(".0", "").replace("_", "-"))
    return {v for v in variants if v}

def score_candidate(fname, protein_name, param_variants):
    """
    Score a candidate filename; higher is better.
    - must start with phi_pairwise_contact_well_ and contain protein_name
    - points if contains param variant
    - extra points if contains decoy/randomization/CPLEX
    - small penalty if contains 'native' only
    """
    score = 0
    if not fname.startswith("phi_pairwise_contact_well_"):
        return -999
    if f"{protein_name}" not in fname:
        return -999

    # param match
    param_matched = any(p in fname for p in param_variants)
    if param_matched:
        score += 50
    # prefer names containing decoy/randomization/CPLEX
    if "randomization" in fname or "CPLEX" in fname or "decoy" in fname or "decoys" in fname:
        score += 30
    # prefer decoys over native
    if "decoy" in fname or "decoys" in fname:
        score += 20
    if "native" in fname:
        score -= 5
    # small tie-breakers
    if fname.endswith("_") or fname.endswith("-"):
        score -= 1
    score += len(fname) % 7  # arbitrary small differentiator
    return score if param_matched else score - 20

def detect_decoy_file(training_set_file, full_name):
    """
    Robustly find the best matching file in ./phis for the given training set and phi params.
    full_name is the 'full_name' returned by get_total_phis_and_parameter_string (may contain
    'phi_pairwise_contact_well-...' format).
    """
    proteins = read_column_from_file(training_set_file, 1)
    if len(proteins) == 0:
        raise RuntimeError("ERROR: training_set_file contains no proteins")

    protein_name = proteins[0].strip()  # e.g., '1urn', '2c4q'
    param_variants = normalize_param_string(full_name)

    # List directory and score candidates
    candidates = []
    try:
        files = os.listdir(phis_directory)
    except Exception as e:
        raise RuntimeError(f"ERROR: cannot list {phis_directory}: {e}")

    for f in files:
        score = score_candidate(f, protein_name, param_variants)
        if score > -100:
            candidates.append((score, f))

    # If no candidates, print debug and raise
    if not candidates:
        print("\nDEBUG: No candidate files found matching basic pattern.")
        print(f"Looking for protein: '{protein_name}' and one of param variants: {sorted(list(param_variants))}")
        print("Files in phis/:")
        for line in sorted(files):
            print("  " + line)
        raise RuntimeError(
            f"ERROR: No decoy file found for protein '{protein_name}' with phi params '{full_name}'"
        )

    # pick the highest scoring candidate
    candidates.sort(reverse=True, key=lambda x: (x[0], x[1]))
    best_score, best_fname = candidates[0]

    # If top candidates are close in score, optionally show warning
    if len(candidates) > 1 and candidates[0][0] == candidates[1][0]:
        print("WARNING: multiple candidate decoy files matched; choosing highest-sorted filename. Candidates:")
        for sc, fn in candidates[:5]:
            print(f"  score={sc:3d}  file={fn}")

    # final check: if best score is low, show debug and fail
    if best_score < 0:
        print("\nDEBUG: Best candidate score is low; showing candidates:")
        for sc, fn in candidates[:10]:
            print(f"  score={sc:3d}  file={fn}")
        raise RuntimeError(
            f"ERROR: No confident decoy file match for protein '{protein_name}' with phi params '{full_name}'"
        )

    return best_fname

# -------------------------------
# Helper functions
# -------------------------------
def sort_eigenvalues_and_eigenvectors(lamb, P):
    idx = np.argsort(lamb)[::-1]  # largest first
    return lamb[idx], P[:, idx]

def estimate_dynamic_cutoff(lamb, half_B, other_half_B, std_half_B, num_decoys):
    cutoffs = []
    total_phis = len(lamb)
    for _ in range(noise_iterations):
        noisy_diag = np.random.normal(
            loc=np.diag(half_B),
            scale=(np.diag(std_half_B) / float(num_decoys) if num_decoys > 0 else np.diag(std_half_B))
        ) - np.diag(other_half_B)

        noisy = np.copy(lamb)
        tail_len = min(len(noisy_diag), len(noisy))
        noisy[-tail_len:] = noisy_diag[-tail_len:]

        dif = np.abs(lamb - noisy) / (np.abs(lamb) + 1e-12)
        idx = np.where(dif > relative_error_threshold)[0]
        cutoffs.append(int(idx[0]) if idx.size else len(lamb)-1)

    return max(3, min(cutoffs))

def calc_A_B_memory_safe(phi_native_mean, phi_decoys):
    total_phis = phi_decoys.shape[1]
    half = np.zeros((total_phis, total_phis), dtype=np.float64)
    std  = np.zeros((total_phis, total_phis), dtype=np.float64)
    oth  = np.zeros((total_phis, total_phis), dtype=np.float64)

    n_rows = phi_decoys.shape[0]
    phis_i = phi_decoys.reshape(n_rows, total_phis, 1)
    phis_j = phi_decoys.reshape(n_rows, 1, total_phis)
    mul = phis_i * phis_j
    half += np.mean(mul, axis=0)
    std  += np.std(mul, axis=0)
    avg = np.mean(phi_decoys, axis=0)
    oth += avg.reshape(total_phis,1) @ avg.reshape(1,total_phis)

    return half, std, oth, avg

def get_filtered_B_inv_lambda_and_P(filtered_lamb, cutoff, P):
    filtered_lamb_f = filtered_lamb.copy()
    filtered_lamb_f[cutoff:] = filtered_lamb_f[cutoff-1]

    with np.errstate(divide='ignore', invalid='ignore'):
        inv_diag = 1.0 / filtered_lamb_f
        inv_diag[~np.isfinite(inv_diag)] = 0.0

    B_inv = P @ np.diag(inv_diag) @ np.linalg.inv(P)
    return B_inv, filtered_lamb_f, P

# -------------------------------
# Core solver
# -------------------------------
def solve(training_set_file, phi_list_file):

    # Read phi list first so we know the parameter string
    phi_list = read_phi_list(phi_list_file)
    training_set = read_column_from_file(training_set_file, 1)
    total_phis, full_name, num_phis = get_total_phis_and_parameter_string(phi_list, training_set)

    # Detect decoy file with both protein + phi params (robust)
    decoy_phi_file = detect_decoy_file(training_set_file, full_name)
    print(f"Using decoy file: {decoy_phi_file}")

    base_name = os.path.splitext(os.path.basename(training_set_file))[0]

    # Native phi summary
    phi_native = np.zeros((len(training_set), total_phis), dtype=np.float64)
    for i, protein in enumerate(training_set):
        phi_native[i] = read_native_phi(protein, phi_list, total_phis)
    phi_native_mean = np.mean(phi_native, axis=0)
    np.savetxt(f"{phis_directory}{base_name}_native_summary.txt", phi_native_mean, fmt='%1.5f')

    # Load decoys
    decoy_path = os.path.join(phis_directory, decoy_phi_file)
    print(f"Loading decoy file path: {decoy_path}")
    phi_decoys = np.loadtxt(decoy_path)
    num_decoys = phi_decoys.shape[0]
    print(f"Detected {num_decoys} decoys in {decoy_phi_file}")
    print(f"phi_native_mean.shape = {phi_native_mean.shape}, phi_decoys.shape = {phi_decoys.shape}")

    # Compute matrices
    half, std, oth, avg_decoy = calc_A_B_memory_safe(phi_native_mean, phi_decoys)
    np.savetxt(f"{phis_directory}{base_name}_decoy_summary.txt", avg_decoy, fmt='%1.5f')

    B = half - oth
    A = avg_decoy - phi_native_mean

    B = (B + B.T) / 2.0
    lamb, P = np.linalg.eig(B)
    lamb, P = sort_eigenvalues_and_eigenvectors(lamb, P)

    cutoff = manual_cutoff if not use_dynamic_cutoff else \
        estimate_dynamic_cutoff(lamb, half, oth, std, num_decoys)
    cutoff = max(1, min(len(lamb)-1, cutoff))

    B_inv, lamb_f, P = get_filtered_B_inv_lambda_and_P(lamb, cutoff, P)
    gamma = B_inv @ A

    out_prefix = f"{gammas_directory}{base_name}_{full_name}"
    np.savetxt(out_prefix + "_gamma", gamma, fmt='%1.5f')
    np.savetxt(out_prefix + "_A", A, fmt='%1.5f')
    np.savetxt(out_prefix + "_B", B, fmt='%1.5f')
    np.savetxt(out_prefix + "_lamb", lamb, fmt='%1.5f')
    np.savetxt(out_prefix + "_lamb_filtered", lamb_f, fmt='%1.5f')
    np.savetxt(out_prefix + "_gamma_filtered", gamma, fmt='%1.5f')

    print(f"\n✨ Gamma computed. Cutoff={cutoff} | Mode={'Dynamic' if use_dynamic_cutoff else 'Manual'}")
    print(f"Saved → {gammas_directory}\n")

# -------------------------------
# Run example
# -------------------------------
if __name__ == "__main__":
    solve("native_trainSetFiles.txt", "phi1_list.txt")
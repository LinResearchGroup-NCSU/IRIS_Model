import numpy as np

gamma_file_name = 'results_phi_gamma/native_trainSetFiles_phi_pairwise_contact_well-9.5_9.5_0.7_10_gamma_filtered'
phi_file_name = 'results_phi_gamma/phi_pairwise_contact_well_native_Rmodified_decoys_CPLEX_randomization_-9.5_9.5_0.7_10'

# Load gamma data
gamma = np.loadtxt(
    gamma_file_name,
    dtype=complex,
    converters={0: lambda s: complex(s.decode().replace('+-', '-'))}
)

# Determine number of phis in first line
with open(phi_file_name, "r") as file:
    first_line = file.readline()
    total_phis = len(first_line.strip().split())

# Count decoys
num_decoys = sum(1 for _ in open(phi_file_name))

# Initialize phi array
phi_i_decoy = np.zeros((num_decoys, total_phis))

# Load phi data
with open(phi_file_name, "r") as file:
    for i_decoy, line in enumerate(file):
        if i_decoy >= num_decoys:
            break
        values = line.strip().split()
        for i_phi, value in enumerate(values):
            phi_i_decoy[i_decoy][i_phi] = float(value)

# Calculate energy for each decoy
e_decoy = np.zeros(num_decoys)
for i_decoy in range(num_decoys):
    e_decoy[i_decoy] = np.dot(gamma, phi_i_decoy[i_decoy])

# Save results
np.savetxt('Energy_mg.txt', e_decoy, fmt='%f', delimiter='\n')
print("Energy calculation complete. Results saved to Energy_mg.txt")

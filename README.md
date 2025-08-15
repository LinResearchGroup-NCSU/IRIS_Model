<h1>
  <img
    src=""
    width="200"
    alt="logo"
    align="middle"
  />
</h1>


The **Interpretable RNA–Protein Interactions Informed by Structure (IRIS) Model** is a computational framework that learns RNA–protein physicochemical interactions by fusing available crystal structures and their associated sequences into an optimized energy model. We show that the model can be used to accurately predict the sequence-specific binding affinities and binding sites of RNA-binding proteins and is transferable across the same protein superfamily. This repository provides a clean implementation of the IRIS model, with training and testing data for selected RNA-binding proteins (PDB IDs: 2c4q). Results are used in **Figure S2** of the IRIS manuscript.

## Features

- Train and visualize IRIS energy models for protein–RNA complexes.
- Generate phi values and calculate binding energies for testing binders.
- Supplementary materials.

## Installation

- **Python 3**: Download from [python.org](https://www.python.org/downloads/).
- **Conda** (Recommended: [Miniconda](https://docs.conda.io/en/latest/miniconda.html))
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   
   ```
1. Clone the repository:
   ```bash
git clone https://github.com/LinResearchGroup-NCSU/IRIS_Model
   ```
2. Create and Activate the Conda Environment:
   ```bash
   cd IRIS_Model
   conda env create -f IRIS.yaml -n IRIS
   conda activate IRIS
   ```
2. Ensure scripts are executable:
   ```bash
   chmod +x *.sh
   ```
  **Note**: Our old version used a single script buildseq.py, from the Modeller package to extract the protein and RNA sequences from the given PDB structure. To avoid further confusion and reduce the environment dependencies required to run the code, we rewrote it using basic Python code. The most time-consuming step of IRIS is generating the phi values for decoy/testing binders. We rewrote the template_evaluate_phi.py by adopting parallel computation from joblib and multiprocessing, which is especially useful when we want to do transcriptome-wide prediction. We also simplified redundant parts of the code and added informative messages during model training and testing to make it more user-friendly.

## Usage

### Training the Energy Model

Train an energy model for three MAX protein complexes (2c4q).

1. **Prepare Input Files**:
   - Create `training/proteinList.txt` with PDB IDs:
     ```
     2c4q
     ```
   - Place PDB files in `training/PDBs`, named `{PDB ID}_modified.pdb, {PDB ID}_Rmodified.pdb,` (e.g., `2c4q_modified.pdb, 2c4q_Rmodified.pdb`). Rename protein chains to **A** and DNA chains to **B** and **C**, respectively.

2. **Run Training**:
   ```bash
   bash train.sh
   ```

3. **Configure Settings**:
   - **Interaction Atoms**: Edit `get_interaction_atom` in `common_function.py` to implement a **coarse-grained interaction scheme**, where RNA bases are represented by the **P** atom (or **C5**), and protein residues are represented by the **Cα (CA)** atom, although **side-chain atoms** may also be used in cases where more detailed interactions are of interest.
   - **Decoy Number**: In `training/optimization/for_bindingE/template/sequences`, the scripts generate_decoy_seq_prot.py and generate_decoy_seq_RNA.py are used to generate protein and RNA decoy sequences for training. The default sizes are 10,000 and 1,000, respectively, which we found to be robust across all tests in our manuscript. Please adjust these values based on your specific needs.
   - **Cutoff Mode**: In `training/optimization/for_training_gamma/optimize_gamma.py`, set cutoff_mode = 25 to retain the first 25 eigenvalues, replacing all others with the 25th eigenvalue. This choice typically depends on the lambda values in `training/optimization/for_training_gamma/gammas/randomized_decoy/native_trainSetFiles_phi_pairwise_contact_well-9.5_9.5_0.7_10_lamb`.

4. **Output**:
   - Results are saved in `training/optimization/for_training_gamma/gammas/randomized_decoy`.
   - Key file: `native_trainSetFiles_phi_pairwise_contact_well-9.5_9.5_0.7_10_gamma_filtered` (filtered energy model).

5. **Visualization**:
   - Run:
     ```bash
     cd training/optimization/for_training_gamma/
     python visualize.py
     ```
   - Plots are saved in `training/optimization/for_training_gamma/visualize`.

### Predicting Protein-RNA Binding Energies

Generate phi values and calculate binding energies for given testing binders (e.g., Max 255 mutated binders testing dataset in Maerkl, S. J et al.).

1. **Prepare Input Files**:
   - Place PDB files in `testing/PDBs`.
   - RNA (testing) sequences are in `testing/sequences`.

2. **Generate Testing Phi**:
   ```bash
   cd testing/
   bash test.sh 2c4q
   ```
   - Output in `testing/phis`:
     - `phi_pairwise_contact_well_native_decoys_CPLEX_randomization_-9.5_9.5_0.7_10` (255 lines for testing sequences).
     - `phi_pairwise_contact_well_native_native_9.5_9.5_0.7_10` (1 line for native structure).

3. **Calculate Binding Energy**:
   - Copy the files to `results_phi_gamma/`:
     ```bash
     cp training/optimization/for_training_gamma/gammas/randomized_decoy/native_trainSetFiles_phi_pairwise_contact_well-9.5_9.5_0.7_10_gamma_filtered results_phi_gamma/
     cp testing/phis/phi_pairwise_contact_well_native_decoys_CPLEX_randomization_-9.5_9.5_0.7_10 results_phi_gamma/
     ```
   
   - Navigate to `IRIS_analysis` and run:
     ```bash
     python energy_calculation.py
     ```
   - Output: `Energy_mg.txt` (predicted energies using **E = γΦ**).

## Supplementary Materials

- **Trained Energy Models**: 
- **Raw Data**: 
- **Processed Published Models**: 

## References

- Zhang, Y., Silvernail, I., Lin, Z., & Lin, X. (2025). Interpretable Protein–DNA Interactions Captured by Structure–Sequence Optimization. eLife, 14, e105565. https://doi.org/10.7554/eLife.105565
- Maerkl, S. J., & Quake, S. R. (2007). A Systems Approach to Measuring the Binding Energy Landscapes of Transcription Factors. *Science*, 315(5809), 233–237. [DOI:10.1126/science.1131007](https://doi.org/10.1126/science.1131007)

## Contact

For questions or support, contact the [Lin Research Group](https://github.com/LinResearchGroup-NCSU) or open an issue on GitHub.

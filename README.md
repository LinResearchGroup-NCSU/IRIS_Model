<p align="center">
  <img src="https://github.com/LinResearchGroup-NCSU/IRIS_Model/blob/main/IRIS_logo.png" alt="IRIS Logo" width="300"/>
</p>

<p align="center"><b>Integrative RNA–Protein Interaction Prediction Informed by Structure and Sequence (IRIS)</b></p>

<p align="center">
  <a href="https://github.com/LinResearchGroup-NCSU/IRIS_Model"><img src="https://img.shields.io/badge/python-3.8-blue.svg" alt="Python 3.8"></a>
  <a href="https://github.com/LinResearchGroup-NCSU/IRIS_Model/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="License"></a>
</p>



## 🚀 Features

  - **Train & Visualize**: Easily train and visualize IRIS energy models for your protein–RNA complexes of interest.
  - **Predict Binding Energies**: Generate feature vectors (Φ values) and calculate binding energies ($E = \\gamma \\Phi$) for novel sequences.
  - **Supplementary Materials**: Access all data and models used in the original manuscript.

-----

## ⚙️ Installation

### Prerequisites

  - **Python 3**: Download from [python.org](https://www.python.org/downloads/).
  - **Conda**: We recommend [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for environment management.

<!-- end list -->

```bash
# Download and install Miniconda (if you don't have it)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### Setup Steps

1.  **Clone the Repository**:

    ```bash
    git clone https://github.com/LinResearchGroup-NCSU/IRIS_Model
    cd IRIS_Model
    ```

2.  **Create and Activate the Conda Environment**:
    This command creates a new environment named `IRIS` with all the required dependencies listed in `IRIS.yaml`.

    ```bash
    conda env create -f IRIS.yaml -n IRIS
    conda activate IRIS
    ```

3.  **Make Scripts Executable**:

    ```bash
    chmod +x *.sh
    ```

-----

## ▶️ Usage

The workflow is divided into two main parts: training an energy model and then using it to predict binding energies.

### 1\. Training the Energy Model

Here, we'll train an energy model for the MS2-RNA complex (PDB ID: `2c4q`) as an example.

#### Step 1: Prepare Input Files

1.  Create a list of PDB IDs for training in `IRIS_Model/training/proteinList.txt`. Each ID should be on a new line.
    ```
    2c4q
    ```
2.  Place your processed PDB files in the `IRIS_Model/training/PDBs/` directory.
      * **Protein File**: `{PDB_ID}_modified.pdb` (e.g., `2c4q_modified.pdb`). **Rename protein chain to A**.
      * **RNA File**: `{PDB_ID}_Rmodified.pdb` (e.g., `2c4q_Rmodified.pdb`). **Rename RNA chains to B and C**.

#### Step 2: Run the Training Script

Execute the main training script. This will perform the optimization and generate the energy model ($\\gamma$ values).

```bash
bash train.sh
```

#### Step 3: Configure Training Settings (Optional)

You can customize the model's behavior by editing the following files:

  * **Interaction Scheme**: In `common_function.py`, modify the `get_interaction_atom` function to define your coarse-graining level.

      * **Default**: RNA bases are represented by the **P** atom, and protein residues by the **Cα (CA)** atom.
      * **Advanced**: You can specify **side-chain atoms** for a more detailed interaction model.

  * **Decoy Sequences**: In the `IRIS_Model/training/optimization/for_bindingE/template/sequences/` directory, the scripts `generate_decoy_seq_prot.py` and `generate_decoy_seq_RNA.py` generate decoy sequences for training.

      * **Default**: 0 protein decoys and 10,000 RNA decoys. These values were found to be robust in our manuscript.

  * **Eigenvalue Cutoff**: In `IRIS_Model/training/optimization/for_training_gamma/optimize_gamma.py`, set `cutoff_mode`.

      * **Example**: `cutoff_mode = 25` retains the top 25 eigenvalues and replaces all others with the 25th eigenvalue. This choice depends on the lambda values found in `.../gammas/randomized_decoy/native_trainSetFiles_..._lamb`.

#### Step 4: Locate and Visualize Output

  * **Primary Output**: The trained energy model is saved at:
    `IRIS_Model/training/optimization/for_training_gamma/gammas/randomized_decoy/native_trainSetFiles_phi_pairwise_contact_well-9.5_9.5_0.7_10_gamma_filtered`

  * **Visualization**: To generate plots of the energy matrix:

    ```bash
    cd IRIS_Model/visualize/
    ```
    There is a jupyter notebook for visualization.

### 2\. Predicting Protein-RNA Binding Energies

Once the model is trained, you can use it to predict binding energies for a set of test sequences.

#### Step 1: Prepare Input Files

1.  Place the testing PDB file(s) in `IRIS_Model/testing/PDBs/`.
2.  List your testing RNA sequences (one per line) in `IRIS_Model/testing/sequences/rna.seq`.

#### Step 2: Generate Feature Vectors (Φ values)

Navigate to the `IRIS_Model/testing/` directory and run the test script with the target PDB ID.

```bash
cd IRIS_Model/testing/
bash test.sh 2c4q
```

This generates the Φ feature vectors in `IRIS_Model/testing/phis/`.

#### Step 3: Calculate Binding Energy

1.  Copy the trained model ($\\gamma$) and the generated feature vectors (Φ) to the analysis directory:

    ```bash
    # Note: Run these commands from the project's root directory
    cp IRIS_Model/training/optimization/for_training_gamma/gammas/randomized_decoy/native_trainSetFiles_phi_pairwise_contact_well-9.5_9.5_0.7_10_gamma_filtered results_phi_gamma/

    cp IRIS_Model/testing/phis/phi_pairwise_contact_well_native_decoys_CPLEX_randomization_-9.5_9.5_0.7_10 results_phi_gamma/
    ```

2.  Run the final calculation script:

    ```bash
    cd IRIS_analysis/
    python energy_calculation.py
    ```

3.  **Final Output**: The predicted binding energies will be in a file named `Energy_mg.txt`, calculated using the equation **E = γΦ**.

-----

## 📚 Supplementary Materials

  - **Trained Energy Models**: Trained models from our manuscript.
  - **Raw Data**: The raw data used for training and validation.
  - **Processed Published Models**: Processed versions of models from related publications.

*(Links to be added upon publication)*

-----

## 📄 Citation

If you use IRIS in your research, please cite our manuscript.
*(Citation details to be added upon publication)*

-----

## 📞 Contact

For questions, bug reports, or support, please open an issue on this GitHub repository or contact the [Lin Research Group](https://github.com/LinResearchGroup-NCSU).

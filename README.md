# tle-disconnectome-reorganization

Analysis code for the study *Individualized Structural Disconnectomes Shape Multimodal Brain Reorganization After Temporal Lobe Epilepsy Surgery*.

## Overview

Temporal lobe epilepsy surgery removes focal tissue and interrupts distributed white matter pathways. This project examines how each patient's estimated structural disconnectome relates to postoperative functional architecture and longitudinal functional and morphometric brain reorganization.

The study includes a primary surgical cohort from Thomas Jefferson University and an independent surgical cohort from Jinling Hospital. Patient-specific white matter disconnection and multimodal imaging measures are represented as interregional relational matrices and compared using representational similarity analysis.

The main analyses test whether disconnectome geometry corresponds to:

1. Postoperative functional connectivity after accounting for preoperative functional connectivity.
2. Interregional dissimilarity in longitudinal functional connectivity change profiles.
3. Interregional differences in longitudinal gray matter volume change.

Spatial sensitivity analyses separately account for Euclidean distance and differences in regional distance from the surgical cavity. Exploratory analyses evaluate associations with seizure outcome, resection volume, and cognitive change.

## Analysis workflow

| Module | Purpose |
| --- | --- |
| `Block1_TLE_surgery_reorganization` | Calculates nodal functional and morphometric reorganization, performs group inference, and generates regional maps. |
| `Block2_WMD_RSA` | Builds patient-level representational dissimilarity matrices and tests correspondence between white matter disconnection and postoperative brain organization. |
| `Block2_WMD_RSA_Sensitivity` | Repeats the representational similarity analyses with separate spatial covariates. |
| `Block3_Clinical_Association` | Evaluates cognitive change and clinical associations with patient-specific white matter disconnection correspondence coefficients. |

Block 1 Step 0 and Step 1 are intentionally excluded from the public repository because they operate on locally prepared restricted imaging derivatives. The public Block 1 workflow therefore begins at Step 2 and requires prepared regional feature inputs.

## Data structure

The analysis expects the following project-relative structure:

```text
Data/
├── atlas/
├── metadata/
├── TJU/
│   ├── nemo/
│   ├── regional_features/
│   └── surgery_mask/
│       ├── native/
│       ├── Orig_1mm/
│       └── xcpd_2mm/
└── JLH/
    ├── nemo/
    ├── regional_features/
    └── surgery_mask/
        ├── native/
        ├── Orig_1mm/
        └── xcpd_2mm/
```

| Directory | Contents |
| --- | --- |
| `Data/atlas` | Parcellation images, lookup tables, regional coordinates, and spatial reference files. |
| `Data/metadata` | Deidentified cohort metadata and analysis variables for TJU patients, TJU healthy participants, and JLH patients. |
| `Data/<Dataset>/nemo` | Patient-specific continuous ChaCo white matter disconnection matrices and associated regional records. |
| `Data/<Dataset>/regional_features` | Derived regional functional and morphometric features together with regional surgical overlap tables. |
| `Data/<Dataset>/surgery_mask/native` | Deidentified surgical cavity masks in their retained source space. |
| `Data/<Dataset>/surgery_mask/Orig_1mm` | Binary surgical cavity masks in the 1 mm analysis space. |
| `Data/<Dataset>/surgery_mask/xcpd_2mm` | Binary surgical cavity masks in the 2 mm functional analysis space. |

The `Data` directory is currently excluded from version control. It is intended for later release after final privacy review, documentation, and confirmation of the applicable sharing terms.

## Repository structure

```text
tle-disconnectome-reorganization/
├── src/
│   ├── pipeline/
│   │   ├── Block1_TLE_surgery_reorganization/
│   │   ├── Block2_WMD_RSA/
│   │   ├── Block2_WMD_RSA_Sensitivity/
│   │   └── Block3_Clinical_Association/
│   └── utils/
├── .gitignore
├── LICENSE
└── README.md
```

## Requirements

This repository does not include an environment lock file or a dependency specification file. The required environment and software must be installed manually as described below.

### MATLAB

The MATLAB analyses require:

1. MATLAB with the Statistics and Machine Learning Toolbox.
2. The Bioinformatics Toolbox for false discovery rate correction with `mafdr`.
3. ENIGMA Toolbox version 2.0.0 for cortical and subcortical visualization.

### Python

The Python scripts require a Conda environment named `tle_rsa`. The project has been checked with the following configuration:

1. Python 3.11.15.
2. NumPy 2.4.6.
3. pandas 3.0.3.
4. SciPy 1.17.1.
5. NiBabel 5.4.2.
6. rsatoolbox version 0.3.2.

rsatoolbox version 0.3.2 is required because the RSA implementation follows its model fitting and prediction behavior. The other listed versions document the checked local configuration rather than a formally locked compatibility range.

### Input preparation software

The prepared imaging derivatives used by the analyses were generated with:

1. fMRIPrep version 25.2.5.
2. XCP-D version 26.1.1.
3. CAT version 26.0.rc4.
4. NeMo version 2.1.

These preprocessing applications are not required to execute the released scripts when the prepared regional and disconnection inputs are already available.

## Running the analyses

Run scripts from the repository root and preserve the numbered order within each module. Required participant-level inputs and atlas resources must be prepared locally before execution.

1. Run the available Block 1 scripts from `Step2` through `Step4` to calculate nodal reorganization, perform group statistics, and generate regional maps from prepared regional features.
2. Run `Block2_WMD_RSA/Step1_Build_RDMs.m` to construct the relational matrices for both cohorts.
3. Run the Block 2 Python scripts for postoperative functional organization and longitudinal change analyses.
4. Run `Block2_WMD_RSA/Step4_Plot_RSA_Results.m` to generate group summaries and figures.
5. Repeat the corresponding sequence in `Block2_WMD_RSA_Sensitivity` for the spatial sensitivity analyses.
6. Run the Block 3 scripts in numerical order for cognitive and clinical analyses in the primary cohort.

The scripts resolve paths relative to the repository root. Generated outputs are written under `results/`, which is excluded from version control.

## Data availability

Participant imaging and clinical data are not distributed with this repository. Availability information for eligible deidentified derivatives will be added when the associated public data record is finalized.

## Citation

If you use this code, please cite:

> Zhang Q, Javidi SS, Cao R, Sperling MR, Zhang Z, Tracy JI. Individualized Structural Disconnectomes Shape Multimodal Brain Reorganization After Temporal Lobe Epilepsy Surgery. Manuscript in preparation.

The citation will be updated when a permanent publication record becomes available.

## License

This project is licensed under the GNU General Public License version 3. See [LICENSE](LICENSE) for details.

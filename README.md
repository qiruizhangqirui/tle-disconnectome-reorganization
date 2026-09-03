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

## Source code guide

### Block 1: postsurgical nodal reorganization

* `src/pipeline/Block1_TLE_surgery_reorganization/Step2_Calculate_Nodal_Reorganization.m`: Calculates patient-level nodal connectivity reconfiguration and signed nodal morphometric reorganization from prepared regional features. It creates Original and Laterality Normalized outputs and standardizes TJU patient values against longitudinal healthy participant change.
* `src/pipeline/Block1_TLE_surgery_reorganization/Step3_Run_Nodal_Reorganization_Statistics.m`: Produces regional summary tables for functional and morphometric reorganization. It analyzes left and right TLE separately in the Original representation, analyzes the pooled cohort after laterality normalization, and applies separate regional FDR correction to the TJU normative results.
* `src/pipeline/Block1_TLE_surgery_reorganization/Step4_Plot_Nodal_Reorganization_Maps.m`: Generates cortical, subcortical, and volumetric maps from the regional results. It also calculates cross-center spatial correlations for laterality-normalized reorganization patterns.

### Block 2: primary white matter disconnection RSA

* `src/pipeline/Block2_WMD_RSA/Step1_Build_RDMs.m`: Constructs five patient-specific representational dissimilarity matrices for preoperative FC, postoperative FC, functional reorganization, morphometric reorganization, and white matter disconnection. Surgery-affected regions and cerebellar regions are excluded before matrices are written.
* `src/pipeline/Block2_WMD_RSA/Step2_Run_PostFC_RSA.py`: Tests the incremental correspondence of white matter disconnection with postoperative FC after accounting for preoperative FC. It estimates whole-brain and network-level coefficients for both cohorts and uses 2,000 Freedman Lane permutations for whole-brain predictor inference.
* `src/pipeline/Block2_WMD_RSA/Step3_Run_Change_RSA.py`: Tests the correspondence of white matter disconnection with functional and morphometric reorganization dissimilarity. It runs whole-brain region-label permutation tests and estimates network-level patient coefficients for both cohorts.
* `src/pipeline/Block2_WMD_RSA/Step4_Plot_RSA_Results.m`: Summarizes patient-level RSA coefficients with center-specific one-sample tests and pooled network estimates. It exports statistical tables, coefficient distributions, unique variance plots, and network forest plots.

### Block 2 sensitivity analysis

* `src/pipeline/Block2_WMD_RSA_Sensitivity/Step1_Build_RDMs.m`: Reconstructs the five primary RDMs and adds interregional Euclidean distance and patient-specific cavity distance difference RDMs. This creates the seven matrices required by the two spatial sensitivity branches.
* `src/pipeline/Block2_WMD_RSA_Sensitivity/Step2_Run_PostFC_RSA.py`: Repeats the postoperative FC RSA with Euclidean distance or cavity distance difference entered as separate covariates. Each branch uses 2,000 Freedman Lane permutations and writes cohort-specific results.
* `src/pipeline/Block2_WMD_RSA_Sensitivity/Step3_Run_Change_RSA.py`: Repeats the functional and morphometric reorganization models with the two spatial covariates evaluated in separate branches. Whole-brain coefficients use 2,000 Freedman Lane permutations.
* `src/pipeline/Block2_WMD_RSA_Sensitivity/Step4_Plot_RSA_Results.m`: Produces center-specific statistical summaries and figures for the Euclidean distance and cavity distance difference sensitivity branches.

### Block 3: cognitive and clinical analyses

* `src/pipeline/Block3_Clinical_Association/Step1_Analyze_Pre_Post_Cognition.m`: Quantifies preoperative to postoperative change across 13 neuropsychological measures using complete paired observations. It applies FDR correction across measures and creates statistical tables and a standardized paired display.
* `src/pipeline/Block3_Clinical_Association/Step2_Run_TJU_Clinical_Associations.m`: Tests 21 associations between three whole-brain WMD correspondence coefficients and seven surgical, clinical, and cognitive variables in the TJU cohort. It applies one FDR correction across all tests and exports tables, individual plots, and a summary heatmap.

### Shared MATLAB utilities

* `src/utils/calculate_box_whisker_limits.m`: Calculates Tukey lower and upper whisker endpoints from finite observations for consistent boxplot rendering.
* `src/utils/calculate_nodal_connectivity_reconfiguration.m`: Calculates nodal connectivity reconfiguration as one minus the correlation between preoperative and postoperative FC profiles after excluding the target ROI and surgery-affected regions.
* `src/utils/calculate_nodal_reorganization_hc_summary.m`: Creates the regional healthy participant summary used for TJU reference mean and standard deviation maps.
* `src/utils/calculate_nodal_reorganization_regional_results.m`: Creates a regional NCR and NMR result table for one fixed cohort and analysis representation. It calculates raw summaries and, when requested, TJU normative inference with separate FDR families.
* `src/utils/calculate_surgery_overlap.m`: Measures the percentage of each atlas ROI intersected by a surgical cavity and classifies partial and threshold-level surgical involvement.
* `src/utils/calculate_welch_two_sample_statistics.m`: Calculates a two-sided Welch test, group difference with confidence interval, and Hedges' g for two independent groups.
* `src/utils/extract_regional_gmv.m`: Integrates modulated CAT gray matter values within each complete atlas ROI and converts the result to cubic millimeters.
* `src/utils/extract_roi_timeseries.m`: Extracts regional mean BOLD time series while applying the fixed spatial support and postoperative cavity exclusion rules.
* `src/utils/generate_group_cavity_mask.m`: Identifies ROIs affected by surgery in at least one patient within a specified cohort, representation, and laterality subgroup.
* `src/utils/load_wmd_csv_matrix.m`: Reads a fixed CocoYeo143 NeMo connectivity CSV and returns its numeric white matter disconnection matrix.
* `src/utils/plot_cortical_overlay.m`: Renders Schaefer cortical values and the group cavity mask as separate surface layers with a shared color scale.
* `src/utils/plot_patient_rdm_montage.m`: Combines all RDM types for one patient into a single diagnostic montage.
* `src/utils/plot_rdm_matrix.m`: Exports one patient RDM as a labeled high-resolution matrix heatmap.
* `src/utils/plot_subcortical_overlay.m`: Renders subcortical values and the group cavity mask as separate mesh layers with a shared color scale.
* `src/utils/plot_wmd_clinical_boxplot.m`: Displays categorical clinical associations with individual WMD coefficients, group distributions, Hedges' g, and the FDR-adjusted P value.
* `src/utils/plot_wmd_clinical_heatmap.m`: Displays the fixed 7 by 3 clinical association effect matrix with FDR-adjusted P values.
* `src/utils/plot_wmd_clinical_scatter.m`: Displays continuous clinical associations with observations, a least-squares line, a confidence band, Spearman rho, and the FDR-adjusted P value.
* `src/utils/read_cocoyeo143_lut.m`: Reads the CocoYeo143 lookup table, preserves atlas order, and assigns cortical network or subcortical labels.
* `src/utils/select_regional_features.m`: Selects the preoperative and postoperative BOLD and GMV arrays for the Original or Flipped LR atlas orientation.
* `src/utils/write_rdm_matrix_csv.m`: Writes a square RDM with explicit ordered ROI identifiers in its first row and column.
* `src/utils/write_surgical_cavity_frequency_outputs.m`: Generates center-level voxelwise cavity frequency maps and regional surgical involvement summaries.

### Shared Python utilities

* `src/utils/wmd_rsa_analysis.py`: Implements RDM standardization, rsatoolbox regression, unique variance calculation, deterministic permutation seeds, whole-brain inference, network-restricted models, and both spatial sensitivity branches.
* `src/utils/wmd_rsa_io.py`: Defines project-relative Block 2 paths, reads cohort identifiers and RDMs, and maps retained ROIs to the eight network labels.
* `src/utils/wmd_rsa_reporting.py`: Enforces the output schemas and writes primary, network-level, and spatial sensitivity RSA result tables.

## Requirements

This repository does not include an environment lock file or a dependency specification file. The required environment and software must be installed manually as described below.

### MATLAB

The MATLAB analyses require:

1. MATLAB with the Statistics and Machine Learning Toolbox.
2. The Bioinformatics Toolbox for false discovery rate correction with `mafdr`.
3. ENIGMA Toolbox version 2.0.0 for cortical and subcortical visualization.

### Python

The Python analyses require:

1. Python 3.11.15.
2. NumPy 2.4.6.
3. pandas 3.0.3.
4. SciPy 1.17.1.
5. NiBabel 5.4.2.
6. rsatoolbox version 0.3.2.

rsatoolbox version 0.3.2 is required because the RSA implementation follows its model fitting and prediction behavior.

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


## Citation

If you use this code, please cite the associated study:

> Qirui Zhang, Sam S Javidi, Ruoyi Cao, Michael R Sperling, Zhiqiang Zhang, and Joseph I Tracy. *Individualized Structural Disconnectomes Shape Multimodal Brain Reorganization After Temporal Lobe Epilepsy Surgery*. Manuscript in preparation.

This section will be updated with the journal citation and DOI when the article is published.

## License

This project is licensed under the GNU General Public License version 3. See [LICENSE](LICENSE) for details.

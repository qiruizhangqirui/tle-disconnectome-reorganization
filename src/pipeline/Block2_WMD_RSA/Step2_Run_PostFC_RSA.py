# %% Purpose and analysis scope
"""Purpose:
Run patient level and network level postoperative FC RSA for TJU and JLH.

The analysis evaluates POST_FC_RDM as a function of PRE_FC_RDM and WMD_RDM.
It estimates the incremental correspondence of WMD after accounting for
preoperative functional organization.

Whole brain predictor significance is evaluated with 2,000 Freedman Lane
reduced model residual permutations. Network models estimate patient level
coefficients for subsequent cohort level inference without first level
permutation testing.

Author: Qirui Zhang
Created: 2026-09-01
"""

# %% Standard library imports
from __future__ import annotations

import sys
from pathlib import Path


# %% Project path and utility imports
PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "src" / "utils"))

from wmd_rsa_analysis import fit_network_post_fc_rsa, fit_post_fc_rsa  # noqa: E402
from wmd_rsa_io import (  # noqa: E402
    BLOCK2_RESULTS,
    NETWORK_NAMES,
    assert_tle_rsa_runtime,
    load_roi_network_map,
    load_subject_ids,
    load_subject_rdm_set,
)
from wmd_rsa_reporting import (  # noqa: E402
    write_network_post_fc_results,
    write_post_fc_results,
)


# %% Fixed analysis parameters
DATASETS = ("TJU", "JLH")
REQUIRED_RDMS = ("POST_FC_RDM", "PRE_FC_RDM", "WMD_RDM")
N_PERMUTATIONS = 2000


# %% Load the approved shared network definition
assert_tle_rsa_runtime()
roi_network_map = load_roi_network_map()


# %% Estimate whole brain and network coefficients for each cohort
patient_rows_by_dataset = {}
network_rows_by_dataset = {}

for dataset in DATASETS:
    subject_ids = load_subject_ids(dataset)
    patient_rows = []
    network_rows = []

    for subject_index, subject_id in enumerate(subject_ids, start=1):
        print(
            f"{dataset} patient {subject_index:02d}/{len(subject_ids):02d}: "
            f"{subject_id}"
        )
        matrices, roi_ids = load_subject_rdm_set(
            dataset,
            subject_id,
            REQUIRED_RDMS,
        )
        patient_rows.append(
            fit_post_fc_rsa(
                dataset,
                subject_id,
                matrices["POST_FC_RDM"],
                matrices["PRE_FC_RDM"],
                matrices["WMD_RDM"],
                roi_ids,
                N_PERMUTATIONS,
            )
        )
        network_rows.extend(
            fit_network_post_fc_rsa(
                dataset,
                subject_id,
                matrices["POST_FC_RDM"],
                matrices["PRE_FC_RDM"],
                matrices["WMD_RDM"],
                roi_ids,
                roi_network_map,
                NETWORK_NAMES,
            )
        )

    patient_rows_by_dataset[dataset] = patient_rows
    network_rows_by_dataset[dataset] = network_rows


# %% Write approved cohort specific result tables
for dataset in DATASETS:
    dataset_result_directory = BLOCK2_RESULTS / dataset
    write_post_fc_results(
        dataset_result_directory / "POST_FC_RSA_Results.csv",
        patient_rows_by_dataset[dataset],
    )
    write_network_post_fc_results(
        dataset_result_directory / "POST_FC_Network_RSA_Results.csv",
        network_rows_by_dataset[dataset],
    )

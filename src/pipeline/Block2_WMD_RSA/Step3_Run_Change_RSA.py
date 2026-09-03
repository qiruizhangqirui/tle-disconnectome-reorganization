# %% Purpose and analysis scope
"""Purpose:
Run patient level and network level longitudinal change RSA for TJU and JLH.

The analysis evaluates two WMD correspondence models:
1. FUNC_REORG_RDM as a function of WMD_RDM.
2. MORPH_REORG_RDM as a function of WMD_RDM.

Whole brain models use 2,000 two sided region label permutations. Network
models estimate patient level coefficients for subsequent cohort level
inference without first level permutation testing.

Author: Qirui Zhang
Created: 2026-09-02
"""

# %% Standard library imports
from __future__ import annotations

import sys
from pathlib import Path


# %% Project path and utility imports
PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "src" / "utils"))

from wmd_rsa_analysis import fit_change_rsa, fit_network_change_rsa  # noqa: E402
from wmd_rsa_io import (  # noqa: E402
    BLOCK2_RESULTS,
    NETWORK_NAMES,
    assert_tle_rsa_runtime,
    load_roi_network_map,
    load_subject_ids,
    load_subject_rdm_set,
)
from wmd_rsa_reporting import (  # noqa: E402
    write_change_results,
    write_network_change_results,
)


# %% Fixed analysis parameters
DATASETS = ("TJU", "JLH")
ANALYSES = (
    ("FUNC_REORG", "FUNC_REORG_RDM"),
    ("MORPH_REORG", "MORPH_REORG_RDM"),
)
REQUIRED_RDMS = tuple(rdm_name for _, rdm_name in ANALYSES) + ("WMD_RDM",)
N_PERMUTATIONS = 2000


# %% Load the approved shared network definition
assert_tle_rsa_runtime()
roi_network_map = load_roi_network_map()


# %% Estimate whole brain and network coefficients for each cohort
patient_rows_by_dataset = {}
network_rows_by_dataset = {}

for dataset in DATASETS:
    subject_ids = load_subject_ids(dataset)
    patient_rows = {analysis: [] for analysis, _ in ANALYSES}
    network_rows = {analysis: [] for analysis, _ in ANALYSES}

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

        for analysis, rdm_name in ANALYSES:
            patient_rows[analysis].append(
                fit_change_rsa(
                    dataset,
                    subject_id,
                    analysis,
                    matrices[rdm_name],
                    matrices["WMD_RDM"],
                    roi_ids,
                    N_PERMUTATIONS,
                )
            )
            network_rows[analysis].extend(
                fit_network_change_rsa(
                    dataset,
                    subject_id,
                    analysis,
                    matrices[rdm_name],
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
    for analysis, _ in ANALYSES:
        write_change_results(
            dataset_result_directory / f"{analysis}_RSA_Results.csv",
            patient_rows_by_dataset[dataset][analysis],
        )
        write_network_change_results(
            dataset_result_directory / f"{analysis}_Network_RSA_Results.csv",
            network_rows_by_dataset[dataset][analysis],
        )

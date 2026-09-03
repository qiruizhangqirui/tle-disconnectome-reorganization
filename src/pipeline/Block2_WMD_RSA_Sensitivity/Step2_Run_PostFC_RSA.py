# %% Purpose and analysis scope
"""Run spatial sensitivity models for postoperative FC RSA.

Each cohort is evaluated with two separate models:
1. POST_FC_RDM ~ PRE_FC_RDM + WMD_RDM + EUCLIDEAN_DISTANCE_RDM
2. POST_FC_RDM ~ PRE_FC_RDM + WMD_RDM + CAVITY_DISTANCE_DIFF_RDM

Whole brain predictor significance is evaluated with 2,000 Freedman Lane
reduced model residual permutations. TJU and JLH are processed sequentially
and written to separate cohort directories.

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

from wmd_rsa_analysis import fit_post_fc_sensitivity_rsa  # noqa: E402
from wmd_rsa_io import (  # noqa: E402
    BLOCK2_SENSITIVITY_RESULTS,
    assert_tle_rsa_runtime,
    load_subject_ids,
    load_subject_rdm_set,
)
from wmd_rsa_reporting import write_post_fc_sensitivity_results  # noqa: E402


# %% Fixed analysis parameters
DATASETS = ("TJU", "JLH")
N_PERMUTATIONS = 2000
BRANCHES = (
    ("Euclidean_Distance", "EUCLIDEAN_DISTANCE_RDM"),
    ("Cavity_Distance_Difference", "CAVITY_DISTANCE_DIFF_RDM"),
)


# %% Validate the approved Python environment
assert_tle_rsa_runtime()


# %% Estimate each sensitivity model in cohort order
rows_by_branch_and_dataset = {}

for branch_name, covariate_name in BRANCHES:
    required_rdms = (
        "POST_FC_RDM",
        "PRE_FC_RDM",
        "WMD_RDM",
        covariate_name,
    )

    for dataset in DATASETS:
        subject_ids = load_subject_ids(dataset)
        intermediate_directory = (
            BLOCK2_SENSITIVITY_RESULTS / dataset / "intermediate"
        )
        rows = []

        for subject_index, subject_id in enumerate(subject_ids, start=1):
            print(
                f"{branch_name} | {dataset} patient "
                f"{subject_index:02d}/{len(subject_ids):02d}: {subject_id}"
            )
            matrices, roi_ids = load_subject_rdm_set(
                dataset,
                subject_id,
                required_rdms,
                intermediate_dir=intermediate_directory,
            )
            rows.append(
                fit_post_fc_sensitivity_rsa(
                    dataset=dataset,
                    subject_id=subject_id,
                    post_rdm=matrices["POST_FC_RDM"],
                    pre_rdm=matrices["PRE_FC_RDM"],
                    wmd_rdm=matrices["WMD_RDM"],
                    cov_rdm=matrices[covariate_name],
                    cov_name=covariate_name,
                    roi_ids=roi_ids,
                    n_permutations=N_PERMUTATIONS,
                )
            )

        rows_by_branch_and_dataset[(branch_name, dataset)] = rows


# %% Write separate branch and cohort result tables
for branch_name, _ in BRANCHES:
    for dataset in DATASETS:
        output_path = (
            BLOCK2_SENSITIVITY_RESULTS
            / branch_name
            / dataset
            / "POST_FC_RSA_Results.csv"
        )
        write_post_fc_sensitivity_results(
            output_path,
            rows_by_branch_and_dataset[(branch_name, dataset)],
        )
        print(f"Saved: {output_path}")

# %% Purpose and analysis scope
"""Purpose:
Run whole brain longitudinal change spatial sensitivity RSA for TJU and JLH.

Each cohort is evaluated with two separate covariate adjusted models for both
functional and morphometric reorganization. One branch adjusts for
EUCLIDEAN_DISTANCE_RDM and the other adjusts for CAVITY_DISTANCE_DIFF_RDM.

Whole brain predictor significance is evaluated with 2,000 Freedman Lane
reduced model residual permutations. Cohorts and sensitivity branches are
written to separate result directories.

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

from wmd_rsa_analysis import fit_change_sensitivity_rsa  # noqa: E402
from wmd_rsa_io import (  # noqa: E402
    BLOCK2_SENSITIVITY_RESULTS,
    assert_tle_rsa_runtime,
    load_subject_ids,
    load_subject_rdm_set,
)
from wmd_rsa_reporting import write_change_sensitivity_results  # noqa: E402


# %% Fixed analysis parameters
DATASETS = ("TJU", "JLH")
BRANCHES = (
    ("Euclidean_Distance", "EUCLIDEAN_DISTANCE_RDM"),
    ("Cavity_Distance_Difference", "CAVITY_DISTANCE_DIFF_RDM"),
)
ANALYSES = (
    ("FUNC_REORG", "FUNC_REORG_RDM"),
    ("MORPH_REORG", "MORPH_REORG_RDM"),
)
N_PERMUTATIONS = 2000


# %% Validate the approved Python environment
assert_tle_rsa_runtime()


# %% Estimate each sensitivity model in cohort order
rows_by_branch_and_dataset = {}

for branch_name, covariate_name in BRANCHES:
    required_rdms = tuple(rdm_name for _, rdm_name in ANALYSES) + (
        "WMD_RDM",
        covariate_name,
    )

    for dataset in DATASETS:
        subject_ids = load_subject_ids(dataset)
        intermediate_directory = (
            BLOCK2_SENSITIVITY_RESULTS / dataset / "intermediate"
        )
        rows = {analysis: [] for analysis, _ in ANALYSES}

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

            for analysis, rdm_name in ANALYSES:
                rows[analysis].append(
                    fit_change_sensitivity_rsa(
                        dataset=dataset,
                        subject_id=subject_id,
                        analysis=analysis,
                        change_rdm=matrices[rdm_name],
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
        output_directory = BLOCK2_SENSITIVITY_RESULTS / branch_name / dataset
        for analysis, _ in ANALYSES:
            output_path = output_directory / f"{analysis}_RSA_Results.csv"
            write_change_sensitivity_results(
                output_path,
                rows_by_branch_and_dataset[(branch_name, dataset)][analysis],
            )
            print(f"Saved: {output_path}")

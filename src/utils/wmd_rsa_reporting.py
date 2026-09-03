"""Result reporting for Block2 WMD RSA."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


POST_FC_COLUMNS = (
    "Dataset",
    "Subject_ID",
    "Analysis",
    "N_ROIs",
    "Full_Model_R2",
    "Standardized_Beta_PRE",
    "Standardized_Beta_WMD",
    "Unique_Delta_R2_PRE",
    "Unique_Delta_R2_WMD",
    "Permutation_P_WMD",
    "Permutation_P_PRE",
)

CHANGE_COLUMNS = (
    "Dataset",
    "Subject_ID",
    "Analysis",
    "N_ROIs",
    "Model_R2",
    "Standardized_Beta_WMD",
    "Permutation_P_WMD",
)

NETWORK_COMMON_COLUMNS = (
    "Dataset",
    "Subject_ID",
    "Analysis",
    "Network",
    "Network_Order",
    "N_ROIs_Total",
    "N_Network_ROIs",
    "N_Anchored_Edges",
    "N_Within_Network_Edges",
    "N_Network_To_Rest_Edges",
)

NETWORK_POST_FC_COLUMNS = NETWORK_COMMON_COLUMNS + (
    "Full_Model_R2",
    "Standardized_Beta_PRE",
    "Standardized_Beta_WMD",
    "Unique_Delta_R2_WMD",
)

NETWORK_CHANGE_COLUMNS = NETWORK_COMMON_COLUMNS + (
    "Model_R2",
    "Standardized_Beta_WMD",
)


def _write(path: Path, rows: list[dict[str, object]], columns: tuple[str, ...]) -> None:
    table = pd.DataFrame(rows)
    table = table.loc[:, columns].sort_values("Subject_ID")
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False, float_format="%.10g")


def _write_network(
    path: Path, rows: list[dict[str, object]], columns: tuple[str, ...]
) -> None:
    table = pd.DataFrame(rows)
    table = table.loc[:, columns].sort_values(["Subject_ID", "Network_Order"])
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False, float_format="%.10g")


def write_post_fc_results(path: Path, rows: list[dict[str, object]]) -> None:
    """Write the patient level POST FC RSA result table."""
    _write(path, rows, POST_FC_COLUMNS)


def write_change_results(path: Path, rows: list[dict[str, object]]) -> None:
    """Write one patient level change RSA result table."""
    _write(path, rows, CHANGE_COLUMNS)


def write_network_post_fc_results(
    path: Path, rows: list[dict[str, object]]
) -> None:
    """Write patient by network POST FC RSA results."""
    _write_network(path, rows, NETWORK_POST_FC_COLUMNS)


def write_network_change_results(
    path: Path, rows: list[dict[str, object]]
) -> None:
    """Write patient by network change RSA results."""
    _write_network(path, rows, NETWORK_CHANGE_COLUMNS)


POST_FC_SENSITIVITY_COLUMNS = (
    "Dataset",
    "Subject_ID",
    "Analysis",
    "N_ROIs",
    "Full_Model_R2",
    "Standardized_Beta_PRE",
    "Standardized_Beta_WMD",
    "Sensitivity_Covariate",
    "Standardized_Beta_Sensitivity_Covariate",
    "Unique_Delta_R2_PRE",
    "Unique_Delta_R2_WMD",
    "Unique_Delta_R2_Sensitivity_Covariate",
    "Permutation_P_PRE",
    "Permutation_P_WMD",
    "Permutation_P_Sensitivity_Covariate",
)

CHANGE_SENSITIVITY_COLUMNS = (
    "Dataset",
    "Subject_ID",
    "Analysis",
    "N_ROIs",
    "Model_R2",
    "Standardized_Beta_WMD",
    "Sensitivity_Covariate",
    "Standardized_Beta_Sensitivity_Covariate",
    "Partial_R2_WMD",
    "Unique_Delta_R2_WMD",
    "Unique_Delta_R2_Sensitivity_Covariate",
    "Permutation_P_WMD",
    "Permutation_P_Sensitivity_Covariate",
)


def write_post_fc_sensitivity_results(
    path: Path, rows: list[dict[str, object]]
) -> None:
    """Write the patient level POST FC sensitivity RSA result table."""
    _write(path, rows, POST_FC_SENSITIVITY_COLUMNS)


def write_change_sensitivity_results(
    path: Path, rows: list[dict[str, object]]
) -> None:
    """Write one patient level change sensitivity RSA result table."""
    _write(path, rows, CHANGE_SENSITIVITY_COLUMNS)

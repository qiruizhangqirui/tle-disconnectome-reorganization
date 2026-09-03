"""Input and path utilities for Block2 WMD RSA."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BLOCK2_RESULTS = PROJECT_ROOT / "results" / "Block2_WMD_RSA"
BLOCK2_SENSITIVITY_RESULTS = PROJECT_ROOT / "results" / "Block2_WMD_RSA_Sensitivity"
NETWORK_NAMES = (
    "Visual",
    "Somatomotor",
    "DorsalAttention",
    "SalienceVentralAttention",
    "Limbic",
    "Control",
    "Default",
    "Subcortical",
)


def assert_tle_rsa_runtime() -> None:
    """Require the approved tle_rsa Python environment."""
    if Path(sys.executable).resolve().parent.name.lower() != "tle_rsa":
        raise RuntimeError("Block2 WMD RSA must run in the tle_rsa environment")


def load_subject_ids(dataset: str) -> list[str]:
    """Load ordered patient identifiers from the fixed metadata sheet."""
    if dataset == "TJU":
        workbook = PROJECT_ROOT / "Data" / "metadata" / "TJU_metadata.xlsx"
        sheet = "TJU_metadata"
    elif dataset == "JLH":
        workbook = PROJECT_ROOT / "Data" / "metadata" / "JLH_metadata.xlsx"
        sheet = "JLH_metadata"
    else:
        raise ValueError(f"Unsupported dataset: {dataset}")
    table = pd.read_excel(workbook, sheet_name=sheet, usecols=["Subject ID"])
    subject_ids = table["Subject ID"].astype(str).str.strip().tolist()
    if len(subject_ids) != len(set(subject_ids)):
        raise ValueError(f"{dataset} metadata contains duplicate Subject ID values")
    return subject_ids


def load_roi_network_map() -> dict[int, str]:
    """Load the fixed eight-network assignment from the atlas lookup table."""
    lut_path = PROJECT_ROOT / "Data" / "atlas" / "cocoyeo143_LUT.txt"
    network_map: dict[int, str] = {}
    cortical_tokens = {
        "Vis": "Visual",
        "SomMot": "Somatomotor",
        "DorsAttn": "DorsalAttention",
        "SalVentAttn": "SalienceVentralAttention",
        "Limbic": "Limbic",
        "Cont": "Control",
        "Default": "Default",
    }
    with lut_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            roi_id = int(fields[0])
            roi_name = fields[1]
            if 1 <= roi_id <= 16:
                network_map[roi_id] = "Subcortical"
            elif ".7Networks_" in roi_name:
                token = roi_name.split(".7Networks_", maxsplit=1)[1].split("_", maxsplit=1)[0]
                network_map[roi_id] = cortical_tokens[token]
    return network_map


def load_rdm_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load one square matrix CSV with explicit row and column ROI IDs."""
    if not path.is_file():
        raise FileNotFoundError(f"Missing required RDM file: {project_relative(path)}")
    table = pd.read_csv(path)
    if table.columns[0] != "ROI_ID":
        raise ValueError(f"First column must be ROI_ID in {project_relative(path)}")
    row_ids = table.iloc[:, 0].to_numpy(dtype=int)
    column_ids = np.asarray([int(value) for value in table.columns[1:]], dtype=int)
    if not np.array_equal(row_ids, column_ids):
        raise ValueError(f"Row and column ROI IDs differ in {project_relative(path)}")
    matrix = table.iloc[:, 1:].to_numpy(dtype=float)
    if matrix.shape != (row_ids.size, row_ids.size):
        raise ValueError(f"Matrix is not square in {project_relative(path)}")
    return matrix, row_ids


def load_subject_rdm_set(
    dataset: str,
    subject_id: str,
    rdm_names: tuple[str, ...],
    intermediate_dir: Path | None = None,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Load one patient's required RDMs and require one identical ROI order.

    If intermediate_dir is provided, it first checks intermediate_dir for each RDM
    (e.g., branch-specific sensitivity covariates). If not found there, it falls back
    to the base Block2_WMD_RSA intermediate directory for baseline RDMs.
    """
    base_directory = BLOCK2_RESULTS / dataset / "intermediate"
    matrices: dict[str, np.ndarray] = {}
    shared_roi_ids: np.ndarray | None = None
    for name in rdm_names:
        filename = f"{subject_id}_{name}.csv"
        if intermediate_dir is not None and (Path(intermediate_dir) / filename).is_file():
            target_path = Path(intermediate_dir) / filename
        elif (base_directory / filename).is_file():
            target_path = base_directory / filename
        else:
            search_paths = [base_directory / filename]
            if intermediate_dir is not None:
                search_paths.insert(0, Path(intermediate_dir) / filename)
            raise FileNotFoundError(f"Missing required RDM file for {subject_id} '{name}'. Checked: {[str(p) for p in search_paths]}")

        matrix, roi_ids = load_rdm_csv(target_path)
        if shared_roi_ids is None:
            shared_roi_ids = roi_ids
        elif not np.array_equal(roi_ids, shared_roi_ids):
            raise ValueError(f"RDM ROI order differs for {dataset} subject {subject_id}")
        matrices[name] = matrix
    if shared_roi_ids is None:
        raise ValueError("rdm_names must not be empty")
    return matrices, shared_roi_ids



def project_relative(path: Path) -> str:
    """Return a POSIX style path relative to the active project root."""
    return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()

"""Core estimation functions for patient level WMD RSA.

The module fits three classes of models:
1. POST FC on PRE FC and WMD, with unique Delta R squared for both predictors.
2. Functional or morphometric change on WMD.
3. The same coefficient models within eight overlapping network anchored masks.

Whole brain coefficients use rsatoolbox 0.3.2 model fitting. Predictor significance
uses ROI label permutations of a complete RDM. POST models use Freedman Lane reduced
model residual permutations, while change models permute the standardized outcome
RDM directly. Network models return patient level coefficients without patient level
permutation inference.

All fitted RDM vectors use population standardization with ``ddof=0``. This matches
the internal correlation scaling used by rsatoolbox 0.3.2. Network edges outside the
active mask are represented as ``NaN`` so that the toolbox fits only shared finite
entries while preserving the full ROI axis.

"""

from __future__ import annotations

import hashlib

import numpy as np
from rsatoolbox.model import ModelWeighted, fit_regress
from rsatoolbox.rdm import RDMs
from rsatoolbox.rdm.rdms import permute_rdms


def _upper_triangle(matrix: np.ndarray) -> np.ndarray:
    """Return the vector representation of one square RDM."""
    matrix_stack = np.asarray(matrix, dtype=float)[None, :, :]
    return RDMs(dissimilarities=matrix_stack).get_vectors()[0]


def _vector_to_matrix(values: np.ndarray) -> np.ndarray:
    """Restore one condensed RDM vector to a symmetric square matrix."""
    return RDMs(dissimilarities=np.asarray(values, dtype=float)).get_matrices()[0]


def _make_rdms(
    matrices: list[np.ndarray], names: list[str], roi_ids: np.ndarray
) -> RDMs:
    """Create an RDMs object with stable matrix and ROI descriptors."""
    return RDMs(
        dissimilarities=np.stack(matrices),
        rdm_descriptors={"RDM_Name": names},
        pattern_descriptors={"ROI_ID": np.asarray(roi_ids, dtype=int)},
    )


def _standardize_matrix(matrix: np.ndarray) -> np.ndarray:
    """Standardize one RDM vector using the rsatoolbox compatible scale."""
    vector = _upper_triangle(matrix)
    standardized = (vector - vector.mean()) / vector.std(ddof=0)
    return _vector_to_matrix(standardized)


def _fit_weights(
    model: ModelWeighted, outcome_matrix: np.ndarray, roi_ids: np.ndarray
) -> np.ndarray:
    """Estimate unnormalized model weights on the correlation scale."""
    outcome = _make_rdms([outcome_matrix], ["OUTCOME"], roi_ids)
    return np.asarray(
        fit_regress(model, outcome, method="corr", normalize=False), dtype=float
    )


def _fit_weighted_model(
    model_name: str,
    predictor_matrices: list[np.ndarray],
    predictor_names: list[str],
    outcome_matrix: np.ndarray,
    roi_ids: np.ndarray,
) -> tuple[ModelWeighted, np.ndarray, np.ndarray]:
    """Build and fit one RSA model and return its fitted vector."""
    model = ModelWeighted(
        model_name,
        _make_rdms(predictor_matrices, predictor_names, roi_ids),
    )
    beta = _fit_weights(model, outcome_matrix, roi_ids)
    return model, beta, model.predict(beta)


def _r_squared(outcome: np.ndarray, fitted: np.ndarray) -> float:
    """Calculate variance explained on the standardized outcome vector."""
    residual_sum_squares = float(np.sum((outcome - fitted) ** 2))
    total_sum_squares = float(np.sum((outcome - outcome.mean()) ** 2))
    return 1.0 - residual_sum_squares / total_sum_squares


def _seed(dataset: str, subject_id: str, analysis: str) -> int:
    """Create a reproducible 32 bit permutation seed for one analysis."""
    label = "|".join((dataset, subject_id, analysis, "WMD")).encode("utf-8")
    return int.from_bytes(hashlib.sha256(label).digest()[:4], byteorder="little")


def _freedman_lane_permutation_p(
    outcome_vector: np.ndarray,
    roi_ids: np.ndarray,
    full_model: ModelWeighted,
    reduced_model: ModelWeighted,
    beta_index: int,
    observed_beta: float,
    n_permutations: int,
    random_seed: int,
) -> float:
    """Test one full model coefficient with Freedman Lane permutations.

    The reduced model is fitted once. Its residual RDM is permuted by applying the
    same random ROI order to rows and columns. Permuted residuals are added to the
    fixed reduced fitted vector before the full model is fitted again. The returned
    value is a two sided Monte Carlo P value with the standard plus one correction.
    """
    outcome_matrix = _vector_to_matrix(outcome_vector)
    reduced_beta = _fit_weights(reduced_model, outcome_matrix, roi_ids)
    reduced_fitted = reduced_model.predict(reduced_beta)
    residual = outcome_vector - reduced_fitted
    residual_rdms = _make_rdms(
        [_vector_to_matrix(residual)], ["REDUCED_RESIDUAL"], roi_ids
    )

    exceedances = 0
    rng = np.random.default_rng(random_seed)
    for _ in range(n_permutations):
        order = rng.permutation(len(roi_ids))
        permuted_residual = permute_rdms(residual_rdms, order).get_vectors()[0]
        permuted_outcome = _vector_to_matrix(reduced_fitted + permuted_residual)
        permuted_beta = _fit_weights(full_model, permuted_outcome, roi_ids)[beta_index]
        if abs(permuted_beta) >= abs(observed_beta):
            exceedances += 1
    return (1.0 + exceedances) / (n_permutations + 1.0)


def _change_wmd_permutation_p(
    change_matrix: np.ndarray,
    roi_ids: np.ndarray,
    model: ModelWeighted,
    observed_beta: float,
    n_permutations: int,
    random_seed: int,
) -> float:
    """Test WMD by permuting the complete standardized change RDM.

    Each random ROI order is applied jointly to rows and columns, which preserves
    the dependency structure induced by shared RDM nodes.
    """
    change_rdms = _make_rdms([change_matrix], ["CHANGE"], roi_ids)
    exceedances = 0
    rng = np.random.default_rng(random_seed)
    for _ in range(n_permutations):
        order = rng.permutation(len(roi_ids))
        permuted_change = permute_rdms(change_rdms, order).get_vectors()[0]
        permuted_change = _vector_to_matrix(permuted_change)
        permuted_beta = _fit_weights(model, permuted_change, roi_ids)[0]
        if abs(permuted_beta) >= abs(observed_beta):
            exceedances += 1
    return (1.0 + exceedances) / (n_permutations + 1.0)


def fit_post_fc_rsa(
    dataset: str,
    subject_id: str,
    post_rdm: np.ndarray,
    pre_rdm: np.ndarray,
    wmd_rdm: np.ndarray,
    roi_ids: np.ndarray,
    n_permutations: int,
) -> dict[str, object]:
    """Fit one whole brain POST FC model and test both predictors."""
    ids = np.asarray(roi_ids, dtype=int).ravel()
    post = _standardize_matrix(post_rdm)
    pre = _standardize_matrix(pre_rdm)
    wmd = _standardize_matrix(wmd_rdm)
    post_vector = _upper_triangle(post)

    # The two single predictor models define the unique Delta R squared terms.
    full_model, full_beta, full_fitted = _fit_weighted_model(
        "POST_FC", [pre, wmd], ["PRE_FC", "WMD"], post, ids
    )
    pre_model, _, pre_fitted = _fit_weighted_model(
        "PRE_FC", [pre], ["PRE_FC"], post, ids
    )
    wmd_model, _, wmd_fitted = _fit_weighted_model(
        "WMD_ONLY", [wmd], ["WMD"], post, ids
    )
    full_r2 = _r_squared(post_vector, full_fitted)
    pre_r2 = _r_squared(post_vector, pre_fitted)
    wmd_r2 = _r_squared(post_vector, wmd_fitted)

    permutation_p_wmd = _freedman_lane_permutation_p(
        post_vector,
        ids,
        full_model,
        pre_model,
        beta_index=1,
        observed_beta=float(full_beta[1]),
        n_permutations=n_permutations,
        random_seed=_seed(dataset, subject_id, "POST_FC"),
    )
    permutation_p_pre = _freedman_lane_permutation_p(
        post_vector,
        ids,
        full_model,
        wmd_model,
        beta_index=0,
        observed_beta=float(full_beta[0]),
        n_permutations=n_permutations,
        random_seed=_seed(dataset, subject_id, "POST_FC_PRE"),
    )
    return {
        "Dataset": dataset,
        "Subject_ID": subject_id,
        "Analysis": "POST_FC",
        "N_ROIs": ids.size,
        "Full_Model_R2": full_r2,
        "Standardized_Beta_PRE": float(full_beta[0]),
        "Standardized_Beta_WMD": float(full_beta[1]),
        "Unique_Delta_R2_PRE": full_r2 - wmd_r2,
        "Unique_Delta_R2_WMD": full_r2 - pre_r2,
        "Permutation_P_WMD": permutation_p_wmd,
        "Permutation_P_PRE": permutation_p_pre,
    }


def fit_change_rsa(
    dataset: str,
    subject_id: str,
    analysis: str,
    change_rdm: np.ndarray,
    wmd_rdm: np.ndarray,
    roi_ids: np.ndarray,
    n_permutations: int,
) -> dict[str, object]:
    """Fit one whole brain change RDM on WMD and permute the outcome RDM."""
    ids = np.asarray(roi_ids, dtype=int).ravel()
    change = _standardize_matrix(change_rdm)
    wmd = _standardize_matrix(wmd_rdm)
    change_vector = _upper_triangle(change)
    model, beta, fitted = _fit_weighted_model(
        analysis, [wmd], ["WMD"], change, ids
    )
    beta_wmd = float(beta[0])
    permutation_p = _change_wmd_permutation_p(
        change,
        ids,
        model,
        observed_beta=beta_wmd,
        n_permutations=n_permutations,
        random_seed=_seed(dataset, subject_id, analysis),
    )

    return {
        "Dataset": dataset,
        "Subject_ID": subject_id,
        "Analysis": analysis,
        "N_ROIs": ids.size,
        "Model_R2": _r_squared(change_vector, fitted),
        "Standardized_Beta_WMD": beta_wmd,
        "Permutation_P_WMD": permutation_p,
    }


def _prepare_masked_rdm(
    all_edges: np.ndarray, mask: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Standardize selected edges and embed them in a masked full RDM."""
    selected = np.asarray(all_edges, dtype=float)[mask]
    standardized = (selected - selected.mean()) / selected.std(ddof=0)
    masked_values = np.full(np.asarray(mask).shape, np.nan, dtype=float)
    masked_values[mask] = standardized
    return standardized, _vector_to_matrix(masked_values)


def _network_edge_sets(
    roi_ids: np.ndarray,
    roi_network_map: dict[int, str],
    network_names: tuple[str, ...],
) -> list[dict[str, object]]:
    """Build overlapping masks containing within network and network to rest edges."""
    ids = np.asarray(roi_ids, dtype=int).ravel()
    upper = np.triu_indices(ids.size, k=1)
    assignments = np.asarray([roi_network_map[int(roi_id)] for roi_id in ids])
    edge_sets = []
    for network_order, network in enumerate(network_names, start=1):
        member = assignments == network
        within = member[upper[0]] & member[upper[1]]
        anchored = member[upper[0]] | member[upper[1]]
        edge_sets.append(
            {
                "Network": network,
                "Network_Order": network_order,
                "Mask": anchored,
                "N_Network_ROIs": int(member.sum()),
                "N_Anchored_Edges": int(anchored.sum()),
                "N_Within_Network_Edges": int(within.sum()),
                "N_Network_To_Rest_Edges": int(anchored.sum() - within.sum()),
            }
        )
    return edge_sets


def fit_network_post_fc_rsa(
    dataset: str,
    subject_id: str,
    post_rdm: np.ndarray,
    pre_rdm: np.ndarray,
    wmd_rdm: np.ndarray,
    roi_ids: np.ndarray,
    roi_network_map: dict[int, str],
    network_names: tuple[str, ...],
) -> list[dict[str, object]]:
    """Estimate POST FC coefficients within each network anchored edge set."""
    ids = np.asarray(roi_ids, dtype=int).ravel()
    post = _standardize_matrix(post_rdm)
    pre = _standardize_matrix(pre_rdm)
    wmd = _standardize_matrix(wmd_rdm)
    post_all = _upper_triangle(post)
    pre_all = _upper_triangle(pre)
    wmd_all = _upper_triangle(wmd)
    edge_sets = _network_edge_sets(ids, roi_network_map, network_names)

    rows: list[dict[str, object]] = []
    for edge_set in edge_sets:
        mask = np.asarray(edge_set["Mask"], dtype=bool)
        # NaNs outside the mask let rsatoolbox retain the complete ROI axis.
        post_vector, post_masked = _prepare_masked_rdm(post_all, mask)
        _, pre_masked = _prepare_masked_rdm(pre_all, mask)
        _, wmd_masked = _prepare_masked_rdm(wmd_all, mask)
        _, full_beta, full_fitted_all = _fit_weighted_model(
            "POST_FC_NETWORK",
            [pre_masked, wmd_masked],
            ["PRE_FC", "WMD"],
            post_masked,
            ids,
        )
        _, _, reduced_fitted_all = _fit_weighted_model(
            "PRE_FC_NETWORK",
            [pre_masked],
            ["PRE_FC"],
            post_masked,
            ids,
        )
        full_fitted = full_fitted_all[mask]
        reduced_fitted = reduced_fitted_all[mask]
        full_r2 = _r_squared(post_vector, full_fitted)
        reduced_r2 = _r_squared(post_vector, reduced_fitted)
        rows.append(
            {
                "Dataset": dataset,
                "Subject_ID": subject_id,
                "Analysis": "POST_FC",
                "Network": edge_set["Network"],
                "Network_Order": edge_set["Network_Order"],
                "N_ROIs_Total": ids.size,
                "N_Network_ROIs": edge_set["N_Network_ROIs"],
                "N_Anchored_Edges": edge_set["N_Anchored_Edges"],
                "N_Within_Network_Edges": edge_set["N_Within_Network_Edges"],
                "N_Network_To_Rest_Edges": edge_set["N_Network_To_Rest_Edges"],
                "Full_Model_R2": full_r2,
                "Standardized_Beta_PRE": float(full_beta[0]),
                "Standardized_Beta_WMD": float(full_beta[1]),
                "Unique_Delta_R2_WMD": full_r2 - reduced_r2,
            }
        )
    return rows


def fit_network_change_rsa(
    dataset: str,
    subject_id: str,
    analysis: str,
    change_rdm: np.ndarray,
    wmd_rdm: np.ndarray,
    roi_ids: np.ndarray,
    roi_network_map: dict[int, str],
    network_names: tuple[str, ...],
) -> list[dict[str, object]]:
    """Estimate one change RSA coefficient within each network anchored edge set."""
    ids = np.asarray(roi_ids, dtype=int).ravel()
    change = _standardize_matrix(change_rdm)
    wmd = _standardize_matrix(wmd_rdm)
    change_all = _upper_triangle(change)
    wmd_all = _upper_triangle(wmd)
    edge_sets = _network_edge_sets(ids, roi_network_map, network_names)

    rows: list[dict[str, object]] = []
    for edge_set in edge_sets:
        mask = np.asarray(edge_set["Mask"], dtype=bool)
        # The same finite edges are supplied to the outcome and predictor RDMs.
        change_vector, change_masked = _prepare_masked_rdm(change_all, mask)
        _, wmd_masked = _prepare_masked_rdm(wmd_all, mask)
        _, beta, fitted_all = _fit_weighted_model(
            f"{analysis}_NETWORK",
            [wmd_masked],
            ["WMD"],
            change_masked,
            ids,
        )
        beta_wmd = float(beta[0])
        fitted = fitted_all[mask]
        rows.append(
            {
                "Dataset": dataset,
                "Subject_ID": subject_id,
                "Analysis": analysis,
                "Network": edge_set["Network"],
                "Network_Order": edge_set["Network_Order"],
                "N_ROIs_Total": ids.size,
                "N_Network_ROIs": edge_set["N_Network_ROIs"],
                "N_Anchored_Edges": edge_set["N_Anchored_Edges"],
                "N_Within_Network_Edges": edge_set["N_Within_Network_Edges"],
                "N_Network_To_Rest_Edges": edge_set["N_Network_To_Rest_Edges"],
                "Model_R2": _r_squared(change_vector, fitted),
                "Standardized_Beta_WMD": beta_wmd,
            }
        )
    return rows


def fit_post_fc_sensitivity_rsa(
    dataset: str,
    subject_id: str,
    post_rdm: np.ndarray,
    pre_rdm: np.ndarray,
    wmd_rdm: np.ndarray,
    cov_rdm: np.ndarray,
    cov_name: str,
    roi_ids: np.ndarray,
    n_permutations: int,
) -> dict[str, object]:
    """Fit whole brain POST FC with one additional spatial covariate.

    Separate reduced models omit PRE FC, WMD, or the sensitivity covariate. They
    provide unique Delta R squared values and the null residuals required to test
    each coefficient in the full model.
    """
    ids = np.asarray(roi_ids, dtype=int).ravel()
    post = _standardize_matrix(post_rdm)
    pre = _standardize_matrix(pre_rdm)
    wmd = _standardize_matrix(wmd_rdm)
    cov = _standardize_matrix(cov_rdm)

    post_vector = _upper_triangle(post)

    full_model, full_beta, full_fitted = _fit_weighted_model(
        "POST_FC",
        [pre, wmd, cov],
        ["PRE_FC", "WMD", cov_name],
        post,
        ids,
    )
    red_wmd_model, _, red_wmd_fitted = _fit_weighted_model(
        "PRE_COV", [pre, cov], ["PRE_FC", cov_name], post, ids
    )
    red_pre_model, _, red_pre_fitted = _fit_weighted_model(
        "WMD_COV", [wmd, cov], ["WMD", cov_name], post, ids
    )
    red_cov_model, _, red_cov_fitted = _fit_weighted_model(
        "PRE_WMD", [pre, wmd], ["PRE_FC", "WMD"], post, ids
    )

    full_r2 = _r_squared(post_vector, full_fitted)
    red_wmd_r2 = _r_squared(post_vector, red_wmd_fitted)
    red_pre_r2 = _r_squared(post_vector, red_pre_fitted)
    red_cov_r2 = _r_squared(post_vector, red_cov_fitted)

    permutation_p_wmd = _freedman_lane_permutation_p(
        post_vector,
        ids,
        full_model,
        red_wmd_model,
        beta_index=1,
        observed_beta=float(full_beta[1]),
        n_permutations=n_permutations,
        random_seed=_seed(dataset, subject_id, f"POST_FC_{cov_name}_WMD"),
    )
    permutation_p_pre = _freedman_lane_permutation_p(
        post_vector,
        ids,
        full_model,
        red_pre_model,
        beta_index=0,
        observed_beta=float(full_beta[0]),
        n_permutations=n_permutations,
        random_seed=_seed(dataset, subject_id, f"POST_FC_{cov_name}_PRE"),
    )
    permutation_p_cov = _freedman_lane_permutation_p(
        post_vector,
        ids,
        full_model,
        red_cov_model,
        beta_index=2,
        observed_beta=float(full_beta[2]),
        n_permutations=n_permutations,
        random_seed=_seed(dataset, subject_id, f"POST_FC_{cov_name}_COV"),
    )

    return {
        "Dataset": dataset,
        "Subject_ID": subject_id,
        "Analysis": "POST_FC",
        "N_ROIs": ids.size,
        "Full_Model_R2": full_r2,
        "Standardized_Beta_PRE": float(full_beta[0]),
        "Standardized_Beta_WMD": float(full_beta[1]),
        "Sensitivity_Covariate": cov_name,
        "Standardized_Beta_Sensitivity_Covariate": float(full_beta[2]),
        "Unique_Delta_R2_PRE": full_r2 - red_pre_r2,
        "Unique_Delta_R2_WMD": full_r2 - red_wmd_r2,
        "Unique_Delta_R2_Sensitivity_Covariate": full_r2 - red_cov_r2,
        "Permutation_P_PRE": permutation_p_pre,
        "Permutation_P_WMD": permutation_p_wmd,
        "Permutation_P_Sensitivity_Covariate": permutation_p_cov,
    }


def fit_change_sensitivity_rsa(
    dataset: str,
    subject_id: str,
    analysis: str,
    change_rdm: np.ndarray,
    wmd_rdm: np.ndarray,
    cov_rdm: np.ndarray,
    cov_name: str,
    roi_ids: np.ndarray,
    n_permutations: int,
) -> dict[str, object]:
    """Fit a whole brain change RDM on WMD and one spatial covariate.

    Predictor specific reduced models provide unique Delta R squared values and
    Freedman Lane residuals for the two coefficient tests.
    """
    ids = np.asarray(roi_ids, dtype=int).ravel()
    change = _standardize_matrix(change_rdm)
    wmd = _standardize_matrix(wmd_rdm)
    cov = _standardize_matrix(cov_rdm)

    change_vector = _upper_triangle(change)

    full_model, full_beta, full_fitted = _fit_weighted_model(
        analysis, [wmd, cov], ["WMD", cov_name], change, ids
    )
    red_wmd_model, _, red_wmd_fitted = _fit_weighted_model(
        "COV_ONLY", [cov], [cov_name], change, ids
    )
    red_cov_model, _, red_cov_fitted = _fit_weighted_model(
        "WMD_ONLY", [wmd], ["WMD"], change, ids
    )

    full_r2 = _r_squared(change_vector, full_fitted)
    red_wmd_r2 = _r_squared(change_vector, red_wmd_fitted)
    red_cov_r2 = _r_squared(change_vector, red_cov_fitted)

    full_sse = float(np.sum((change_vector - full_fitted) ** 2))
    red_wmd_sse = float(np.sum((change_vector - red_wmd_fitted) ** 2))

    permutation_p_wmd = _freedman_lane_permutation_p(
        change_vector,
        ids,
        full_model,
        red_wmd_model,
        beta_index=0,
        observed_beta=float(full_beta[0]),
        n_permutations=n_permutations,
        random_seed=_seed(dataset, subject_id, f"{analysis}_{cov_name}_WMD"),
    )
    permutation_p_cov = _freedman_lane_permutation_p(
        change_vector,
        ids,
        full_model,
        red_cov_model,
        beta_index=1,
        observed_beta=float(full_beta[1]),
        n_permutations=n_permutations,
        random_seed=_seed(dataset, subject_id, f"{analysis}_{cov_name}_COV"),
    )

    return {
        "Dataset": dataset,
        "Subject_ID": subject_id,
        "Analysis": analysis,
        "N_ROIs": ids.size,
        "Model_R2": full_r2,
        "Standardized_Beta_WMD": float(full_beta[0]),
        "Sensitivity_Covariate": cov_name,
        "Standardized_Beta_Sensitivity_Covariate": float(full_beta[1]),
        "Partial_R2_WMD": max(0.0, (red_wmd_sse - full_sse) / red_wmd_sse),
        "Unique_Delta_R2_WMD": full_r2 - red_wmd_r2,
        "Unique_Delta_R2_Sensitivity_Covariate": full_r2 - red_cov_r2,
        "Permutation_P_WMD": permutation_p_wmd,
        "Permutation_P_Sensitivity_Covariate": permutation_p_cov,
    }

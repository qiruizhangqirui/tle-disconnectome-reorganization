function regional_results = calculate_nodal_reorganization_regional_results( ...
    ncr_table, nmr_table, group_mask_table, dataset, analysis_mode, ...
    tle_laterality, minimum_patient_count, alpha, run_inference)
%CALCULATE_NODAL_REORGANIZATION_REGIONAL_RESULTS Summarizes regional NCR and NMR.
%
% Purpose:
% This utility creates one wide regional result table for a fixed TLE group.
% Raw summaries use all finite nonaffected patient values. TJU normative
% inference uses the subset with finite Normative Z values and positive
% healthy control reference standard deviations. Benjamini Hochberg FDR is
% applied separately to NCR and NMR across eligible nonmasked ROIs.
%
% Author: Qirui Zhang
% Created: 1 September 2026

arguments
    ncr_table table
    nmr_table table
    group_mask_table table
    dataset (1, 1) string
    analysis_mode (1, 1) string
    tle_laterality (1, 1) string
    minimum_patient_count (1, 1) double
    alpha (1, 1) double
    run_inference (1, 1) logical
end

metric_tables = {ncr_table; nmr_table};
raw_value_columns = [ ...
    "Nodal_Connectivity_Reconfiguration"; ...
    "Nodal_Morphometric_Reorganization_Percent"];
reference_sd_columns = [ ...
    "HC_Reference_SD"; ...
    "HC_Reference_SD_Percent"];

patient_rows = string(ncr_table.Group) == "TLE";
if tle_laterality ~= "ALL"
    patient_rows = patient_rows & ...
        upper(string(ncr_table.TLE_Laterality)) == tle_laterality;
end

subject_ids = unique(string(ncr_table.Subject_ID(patient_rows)), 'stable');
first_subject_id = subject_ids(1);
roi_reference = ncr_table( ...
    patient_rows & string(ncr_table.Subject_ID) == first_subject_id, ...
    {'ROI_ID', 'ROI_Name', 'Network'});

roi_ids = roi_reference.ROI_ID;
roi_names = string(roi_reference.ROI_Name);
roi_networks = string(roi_reference.Network);
n_rois = numel(roi_ids);
n_metrics = numel(metric_tables);

n_raw_valid = zeros(n_rois, n_metrics);
mean_raw = nan(n_rois, n_metrics);
sd_raw = nan(n_rois, n_metrics);
n_normative_valid = zeros(n_rois, n_metrics);
mean_normative_z = nan(n_rois, n_metrics);
sd_normative_z = nan(n_rois, n_metrics);
one_sample_t = nan(n_rois, n_metrics);
raw_p = nan(n_rois, n_metrics);
fdr_p = nan(n_rois, n_metrics);
group_masked = false(n_rois, 1);

% Raw values and normative values use separate validity masks so that the
% healthy control reference does not alter the raw regional mean.
for roi_index = 1:n_rois
    mask_row = group_mask_table.ROI_ID == roi_ids(roi_index);
    group_masked(roi_index) = group_mask_table.Group_Masked(mask_row);

    for metric_index = 1:n_metrics
        current_table = metric_tables{metric_index};
        current_group_rows = string(current_table.Group) == "TLE";
        if tle_laterality ~= "ALL"
            current_group_rows = current_group_rows & ...
                upper(string(current_table.TLE_Laterality)) == tle_laterality;
        end

        current_roi_rows = current_group_rows & ...
            current_table.ROI_ID == roi_ids(roi_index);
        current_raw = current_table.(char(raw_value_columns(metric_index)));
        current_surgery_affected = logical(current_table.Surgery_Affected);
        raw_valid_rows = current_roi_rows & ...
            ~current_surgery_affected & isfinite(current_raw);
        raw_values = current_raw(raw_valid_rows);

        n_raw_valid(roi_index, metric_index) = numel(raw_values);
        if ~group_masked(roi_index) && ~isempty(raw_values)
            mean_raw(roi_index, metric_index) = mean(raw_values);
            sd_raw(roi_index, metric_index) = std(raw_values, 0);
        end

        if ~run_inference
            continue;
        end

        current_normative_z = current_table.Normative_Z;
        current_reference_sd = current_table.( ...
            char(reference_sd_columns(metric_index)));
        normative_valid_rows = raw_valid_rows & ...
            isfinite(current_normative_z) & ...
            isfinite(current_reference_sd) & current_reference_sd > 0;
        normative_values = current_normative_z(normative_valid_rows);

        n_normative_valid(roi_index, metric_index) = ...
            numel(normative_values);
        if ~group_masked(roi_index) && ~isempty(normative_values)
            mean_normative_z(roi_index, metric_index) = ...
                mean(normative_values);
            sd_normative_z(roi_index, metric_index) = ...
                std(normative_values, 0);
        end

        if group_masked(roi_index) || ...
                numel(normative_values) < minimum_patient_count
            continue;
        end

        % A one sample t test is undefined when every Normative Z value is
        % identical, so the inferential fields remain NaN for that ROI.
        if std(normative_values, 0) <= eps(max(abs(normative_values)))
            continue;
        end

        [~, current_p, ~, current_statistics] = ttest( ...
            normative_values, 0, 'Alpha', alpha, 'Tail', 'both');
        one_sample_t(roi_index, metric_index) = current_statistics.tstat;
        raw_p(roi_index, metric_index) = current_p;
    end
end

% Each metric defines an independent FDR family containing only eligible
% nonmasked ROIs with a finite one sample test P value.
if run_inference
    for metric_index = 1:n_metrics
        family_rows = ~group_masked & isfinite(raw_p(:, metric_index));
        if any(family_rows)
            fdr_p(family_rows, metric_index) = mafdr( ...
                raw_p(family_rows, metric_index), 'BHFDR', true);
        end
    end
end

n_subjects_total = repmat(numel(subject_ids), n_rois, 1);
regional_results = table( ...
    repmat(dataset, n_rois, 1), ...
    repmat(analysis_mode, n_rois, 1), ...
    repmat("TLE", n_rois, 1), ...
    repmat(tle_laterality, n_rois, 1), ...
    roi_ids, roi_names, roi_networks, group_masked, n_subjects_total, ...
    n_raw_valid(:, 1), mean_raw(:, 1), sd_raw(:, 1), ...
    n_raw_valid(:, 2), mean_raw(:, 2), sd_raw(:, 2), ...
    'VariableNames', { ...
    'Dataset', 'Analysis_Mode', 'Group', 'TLE_Laterality', ...
    'ROI_ID', 'ROI_Name', 'Network', 'Group_Masked', ...
    'N_Subjects_Total', ...
    'NCR_N_Raw_Valid', 'NCR_Mean_Raw', 'NCR_SD_Raw', ...
    'NMR_N_Raw_Valid', 'NMR_Mean_Raw_Percent', ...
    'NMR_SD_Raw_Percent'});

if run_inference
    regional_results = addvars( ...
        regional_results, ...
        n_normative_valid(:, 1), mean_normative_z(:, 1), ...
        sd_normative_z(:, 1), one_sample_t(:, 1), raw_p(:, 1), fdr_p(:, 1), ...
        n_normative_valid(:, 2), mean_normative_z(:, 2), ...
        sd_normative_z(:, 2), one_sample_t(:, 2), raw_p(:, 2), fdr_p(:, 2), ...
        'NewVariableNames', { ...
        'NCR_N_Normative_Valid', 'NCR_Mean_Normative_Z', ...
        'NCR_SD_Normative_Z', 'NCR_One_Sample_T', ...
        'NCR_Raw_P', 'NCR_FDR_P', ...
        'NMR_N_Normative_Valid', 'NMR_Mean_Normative_Z', ...
        'NMR_SD_Normative_Z', 'NMR_One_Sample_T', ...
        'NMR_Raw_P', 'NMR_FDR_P'});
end

end

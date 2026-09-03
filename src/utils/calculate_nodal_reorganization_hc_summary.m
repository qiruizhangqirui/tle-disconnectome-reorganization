function hc_summary = calculate_nodal_reorganization_hc_summary( ...
    ncr_table, nmr_table)
%CALCULATE_NODAL_REORGANIZATION_HC_SUMMARY Summarizes TJU HC raw changes.
%
% Purpose:
% This utility creates the wide regional healthy control summary used for
% TJU Original mean and standard deviation maps. Each row represents one
% retained ROI and contains separate NCR and NMR descriptive statistics.
%
% Author: Qirui Zhang
% Created: 1 September 2026

arguments
    ncr_table table
    nmr_table table
end

hc_rows = string(ncr_table.Group) == "HC";
first_hc_id = string(ncr_table.Subject_ID(find(hc_rows, 1)));
roi_reference = ncr_table( ...
    hc_rows & string(ncr_table.Subject_ID) == first_hc_id, ...
    {'ROI_ID', 'ROI_Name', 'Network'});

roi_ids = roi_reference.ROI_ID;
roi_names = string(roi_reference.ROI_Name);
roi_networks = string(roi_reference.Network);
n_rois = numel(roi_ids);
n_subjects_total = numel(unique(string(ncr_table.Subject_ID(hc_rows))));

ncr_n_valid = zeros(n_rois, 1);
ncr_mean = nan(n_rois, 1);
ncr_sd = nan(n_rois, 1);
nmr_n_valid = zeros(n_rois, 1);
nmr_mean = nan(n_rois, 1);
nmr_sd = nan(n_rois, 1);

% HC summaries use finite raw values because Normative Z scores are defined
% only for patients relative to the healthy control reference distribution.
for roi_index = 1:n_rois
    ncr_roi_rows = hc_rows & ncr_table.ROI_ID == roi_ids(roi_index);
    ncr_values = ncr_table.Nodal_Connectivity_Reconfiguration(ncr_roi_rows);
    ncr_values = ncr_values(isfinite(ncr_values));
    ncr_n_valid(roi_index) = numel(ncr_values);
    ncr_mean(roi_index) = mean(ncr_values);
    ncr_sd(roi_index) = std(ncr_values, 0);

    nmr_roi_rows = string(nmr_table.Group) == "HC" & ...
        nmr_table.ROI_ID == roi_ids(roi_index);
    nmr_values = ...
        nmr_table.Nodal_Morphometric_Reorganization_Percent(nmr_roi_rows);
    nmr_values = nmr_values(isfinite(nmr_values));
    nmr_n_valid(roi_index) = numel(nmr_values);
    nmr_mean(roi_index) = mean(nmr_values);
    nmr_sd(roi_index) = std(nmr_values, 0);
end

hc_summary = table( ...
    repmat("TJU", n_rois, 1), ...
    repmat("Original", n_rois, 1), ...
    repmat("HC", n_rois, 1), ...
    roi_ids, roi_names, roi_networks, ...
    repmat(n_subjects_total, n_rois, 1), ...
    ncr_n_valid, ncr_mean, ncr_sd, ...
    nmr_n_valid, nmr_mean, nmr_sd, ...
    'VariableNames', { ...
    'Dataset', 'Analysis_Mode', 'Group', ...
    'ROI_ID', 'ROI_Name', 'Network', 'N_Subjects_Total', ...
    'NCR_N_Raw_Valid', 'NCR_Mean_Raw', 'NCR_SD_Raw', ...
    'NMR_N_Raw_Valid', 'NMR_Mean_Raw_Percent', ...
    'NMR_SD_Raw_Percent'});

end

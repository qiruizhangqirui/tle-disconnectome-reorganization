function group_mask_table = generate_group_cavity_mask( ...
    overlap_table, dataset, analysis_mode, subgroup_laterality)
%GENERATE_GROUP_CAVITY_MASK Computes group-level surgery cavity mask.
%
% Purpose:
% This utility identifies every ROI affected by surgery in at least one
% patient within a fixed dataset, analysis mode, and laterality subgroup.
% The input Surgery_ROI_Overlap table follows the approved fixed schema.
%
% Syntax:
%   group_mask_table = generate_group_cavity_mask( ...
%       overlap_table, dataset, analysis_mode, subgroup_laterality)
%
% Inputs:
%   overlap_table       - Table containing subject-level surgery overlap
%                         (e.g., Surgery_ROI_Overlap.csv).
%   dataset             - String ("TJU" or "JLH").
%   analysis_mode       - String ("Original" or "Laterality_Normalized").
%   subgroup_laterality - Optional string ("LEFT", "RIGHT", or "ALL"). Default "ALL".
%
% Outputs:
%   group_mask_table    - Table with columns:
%                         Dataset, Analysis_Mode, Subgroup, ROI_ID, ROI_Name,
%                         N_Patients_Total, N_Patients_Surgery_Affected,
%                         Patient_Involvement_Percent, Group_Masked.
%
% Author: Qirui Zhang
% Created: 1 September 2026

arguments
    overlap_table table
    dataset (1, 1) string
    analysis_mode (1, 1) string
    subgroup_laterality (1, 1) string = "ALL"
end

match_rows = string(overlap_table.Dataset) == dataset & ...
    string(overlap_table.Analysis_Mode) == analysis_mode;

if subgroup_laterality ~= "ALL"
    match_rows = match_rows & ...
        upper(string(overlap_table.TLE_Laterality)) == upper(subgroup_laterality);
end

filtered_table = overlap_table(match_rows, :);
patient_ids = unique(string(filtered_table.Subject_ID), 'stable');
n_patients_total = numel(patient_ids);

if n_patients_total == 0
    error( ...
        'Block1:GroupMaskNoPatients', ...
        'No patients found matching dataset=%s, mode=%s, laterality=%s.', ...
        dataset, analysis_mode, subgroup_laterality);
end

roi_reference = unique( ...
    filtered_table(:, {'ROI_ID', 'ROI_Name'}), ...
    'rows', 'stable');
n_rois = height(roi_reference);

group_dataset = repmat(dataset, n_rois, 1);
group_mode = repmat(analysis_mode, n_rois, 1);
group_subgroup = repmat(subgroup_laterality, n_rois, 1);
roi_ids = roi_reference.ROI_ID;
roi_names = string(roi_reference.ROI_Name);
n_affected = zeros(n_rois, 1);
involvement_pct = zeros(n_rois, 1);
group_masked = false(n_rois, 1);

for roi_index = 1:n_rois
    current_roi_id = roi_ids(roi_index);
    current_roi_rows = filtered_table.ROI_ID == current_roi_id;
    current_affected_count = sum(logical(filtered_table.Surgery_Affected(current_roi_rows)));
    n_affected(roi_index) = current_affected_count;
    involvement_pct(roi_index) = (double(current_affected_count) / double(n_patients_total)) * 100;
    group_masked(roi_index) = current_affected_count >= 1;
end

group_mask_table = table( ...
    group_dataset, group_mode, group_subgroup, ...
    roi_ids, roi_names, ...
    repmat(n_patients_total, n_rois, 1), ...
    n_affected, involvement_pct, group_masked, ...
    'VariableNames', { ...
    'Dataset', 'Analysis_Mode', 'Subgroup', ...
    'ROI_ID', 'ROI_Name', ...
    'N_Patients_Total', 'N_Patients_Surgery_Affected', ...
    'Patient_Involvement_Percent', 'Group_Masked'});

end

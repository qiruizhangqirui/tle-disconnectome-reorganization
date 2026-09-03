function [fc_profile_correlation_r, ncr, n_profile_edges] = ...
    calculate_nodal_connectivity_reconfiguration( ...
    pre_fc, post_fc, roi_ids, surgery_affected_roi_ids)
% Purpose:
% Calculate nodal connectivity reconfiguration as one minus the correlation
% between preoperative and postoperative retained-region FC profiles.
%
% Author: Qirui Zhang
% Created: 31 August 2026
%
%CALCULATE_NODAL_CONNECTIVITY_RECONFIGURATION Calculate nodal NCR as 1-r.
%
% Surgery-affected target ROIs are not estimated. For each retained target
% ROI, the PRE and POST profiles exclude the target itself and every
% surgery-affected ROI.
%
% pre_fc and post_fc are ordered 116 by 116 Pearson FC matrices. roi_ids is
% the corresponding 116 by 1 ROI identifier vector.
% surgery_affected_roi_ids identifies target ROIs and profile edges excluded
% for the current patient.
%
% fc_profile_correlation_r contains the PRE to POST Pearson profile
% correlation for each retained target ROI. ncr equals one minus this
% correlation, with higher values indicating greater profile reconfiguration.
% n_profile_edges records the number of retained edges in each correlation.
% Surgery affected targets remain NaN with zero profile edges. Nonfinite
% profile values are not omitted.
%
% This function calculates a subject level descriptive measure and performs
% no hypothesis test or multiple testing correction.

arguments
    pre_fc (:, :) double
    post_fc (:, :) double
    roi_ids (:, 1) double
    surgery_affected_roi_ids (:, 1) double
end

% Initialize ordered output vectors and the retained ROI mask.
n_rois = numel(roi_ids);
affected_mask = ismember(roi_ids, surgery_affected_roi_ids);
retained_mask = ~affected_mask;
fc_profile_correlation_r = nan(n_rois, 1);
ncr = nan(n_rois, 1);
n_profile_edges = zeros(n_rois, 1);

% Calculate retained edge profile correlations for each target ROI.
for roi_index = 1:n_rois
    if affected_mask(roi_index)
        continue;
    end

    profile_mask = retained_mask;
    profile_mask(roi_index) = false;
    n_profile_edges(roi_index) = sum(profile_mask);

    pre_profile = pre_fc(roi_index, profile_mask)';
    post_profile = post_fc(roi_index, profile_mask)';
    correlation_matrix = corrcoef(pre_profile, post_profile);
    fc_profile_correlation_r(roi_index) = correlation_matrix(1, 2);
    ncr(roi_index) = 1 - fc_profile_correlation_r(roi_index);
end
end

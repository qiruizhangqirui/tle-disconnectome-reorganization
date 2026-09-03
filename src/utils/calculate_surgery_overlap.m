function [overlap_percent, roi_voxel_counts, intersection_voxel_counts, ...
    surgery_affected, partially_affected] = calculate_surgery_overlap( ...
    surgery_mask_path, atlas_data, roi_ids, threshold_percent)
% Purpose:
% Calculate regional surgical cavity overlap and classify partially affected
% and surgery-affected regions using the approved percentage threshold.
%
% Author: Qirui Zhang
% Created: 31 August 2026
%
%CALCULATE_SURGERY_OVERLAP Calculate surgery overlap for each atlas ROI.
%
% The input mask is thresholded at 0.5 on the atlas grid. A region is
% affected when overlap meets or exceeds the regional threshold. A region
% is partially affected when overlap is greater than zero but remains below
% the regional threshold.

arguments
    surgery_mask_path (1, 1) string
    atlas_data
    roi_ids (:, 1) double
    threshold_percent (1, 1) double {mustBeNonnegative}
end

mask_data = double(niftiread(surgery_mask_path));
mask_data = squeeze(mask_data);
mask_data = mask_data >= 0.5;
if ~isequal(size(mask_data), size(atlas_data))
    error('Block1:SurgeryOverlapMaskSizeMismatch', ...
        'The surgery overlap mask must match the surgery atlas size: %s', ...
        surgery_mask_path);
end

n_rois = numel(roi_ids);
roi_voxel_counts = zeros(n_rois, 1);
intersection_voxel_counts = zeros(n_rois, 1);
overlap_percent = zeros(n_rois, 1);

for roi_index = 1:n_rois
    roi_mask = atlas_data == roi_ids(roi_index);
    roi_voxel_counts(roi_index) = sum(roi_mask(:));
    intersection_voxel_counts(roi_index) = ...
        sum(mask_data(roi_mask));
    overlap_percent(roi_index) = ...
        100 * intersection_voxel_counts(roi_index) ...
        / roi_voxel_counts(roi_index);
end

surgery_affected = overlap_percent >= threshold_percent;
partially_affected = ...
    overlap_percent > 0 & overlap_percent < threshold_percent;
end

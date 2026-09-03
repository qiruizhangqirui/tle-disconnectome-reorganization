function [time_series, total_roi_voxel_counts, used_roi_voxel_counts, ...
    excluded_roi_voxel_counts] = extract_roi_timeseries( ...
    bold_path, atlas_data, atlas_info, roi_ids, ...
    voxel_exclusion_mask, expected_timepoints)
% Purpose:
% Extract regional mean BOLD time series from a fixed atlas while applying
% the approved spatial support and postoperative cavity exclusion rules.
%
% Author: Qirui Zhang
% Created: 31 August 2026
%
%EXTRACT_ROI_TIMESERIES Extract mean regional signals from a 4D image.
%
% Each session is extracted independently. The spatial support is fixed from
% the first BOLD frame. Atlas voxels excluded by the surgery mask or equal to
% zero or NaN in the first frame are removed. Later zero values are retained,
% while later NaN values are omitted from the regional mean at the affected
% time point. No time points are removed.

arguments
    bold_path (1, 1) string
    atlas_data
    atlas_info
    roi_ids (:, 1) double
    voxel_exclusion_mask
    expected_timepoints (1, 1) double {mustBePositive, mustBeInteger}
end

bold_info = niftiinfo(bold_path);
bold_data = double(niftiread(bold_info));

if ~isequal(size(bold_data, 1), size(atlas_data, 1)) || ...
        ~isequal(size(bold_data, 2), size(atlas_data, 2)) || ...
        ~isequal(size(bold_data, 3), size(atlas_data, 3))
    error('Block1:BOLDSizeMismatch', ...
        'BOLD image size does not match the atlas: %s', bold_path);
end

grid_tolerance = 1e-5;
if ~isequal(double(bold_info.ImageSize(1:3)), ...
        double(atlas_info.ImageSize(1:3)))
    error('Block1:BOLDGridImageSizeMismatch', ...
        'BOLD image size does not match the atlas grid: %s', bold_path);
end
if any(abs(double(bold_info.PixelDimensions(1:3)) ...
        - double(atlas_info.PixelDimensions(1:3))) > grid_tolerance)
    error('Block1:BOLDGridVoxelSizeMismatch', ...
        'BOLD voxel size does not match the atlas grid: %s', bold_path);
end
if any(abs(double(bold_info.Transform.T(:)) ...
        - double(atlas_info.Transform.T(:))) > grid_tolerance)
    error('Block1:BOLDGridAffineMismatch', ...
        'BOLD affine does not match the atlas grid: %s', bold_path);
end

n_timepoints = size(bold_data, 4);
if n_timepoints ~= expected_timepoints
    error('Block1:BOLDTimepointCountMismatch', ...
        'Expected %d BOLD timepoints but found %d: %s', ...
        expected_timepoints, n_timepoints, bold_path);
end
n_rois = numel(roi_ids);
time_series = nan(n_timepoints, n_rois);
total_roi_voxel_counts = zeros(n_rois, 1);
used_roi_voxel_counts = zeros(n_rois, 1);
excluded_roi_voxel_counts = zeros(n_rois, 1);
voxel_by_time = reshape(bold_data, [], n_timepoints);
atlas_vector = atlas_data(:);

if isempty(voxel_exclusion_mask)
    voxel_exclusion_mask = false(size(atlas_data));
end
exclusion_vector = voxel_exclusion_mask(:);

for roi_index = 1:n_rois
    total_roi_mask = atlas_vector == roi_ids(roi_index);
    total_roi_voxel_counts(roi_index) = sum(total_roi_mask);
    candidate_roi_mask = total_roi_mask & ~exclusion_vector;
    candidate_roi_values = voxel_by_time(candidate_roi_mask, :);
    first_frame_values = candidate_roi_values(:, 1);
    valid_signal_voxels = first_frame_values ~= 0 ...
        & ~isnan(first_frame_values);
    roi_values = candidate_roi_values(valid_signal_voxels, :);
    used_roi_voxel_counts(roi_index) = sum(valid_signal_voxels);
    excluded_roi_voxel_counts(roi_index) = ...
        total_roi_voxel_counts(roi_index) ...
        - used_roi_voxel_counts(roi_index);

    if isempty(roi_values)
        error('Block1:NoValidROIVoxels', ...
            'ROI %d has no nonzero non-NaN voxel in the first frame: %s', ...
            roi_ids(roi_index), bold_path);
    end

    roi_signal = mean(roi_values, 1, 'omitnan')';
    if any(~isfinite(roi_signal))
        error('Block1:NonfiniteROISignal', ...
            'ROI %d has a nonfinite regional BOLD value: %s', ...
            roi_ids(roi_index), bold_path);
    end
    time_series(:, roi_index) = roi_signal;
end
end

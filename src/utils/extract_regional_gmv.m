function [gmv_mm3, total_roi_voxel_counts, used_roi_voxel_counts, ...
    excluded_roi_voxel_counts] = extract_regional_gmv( ...
    gmv_path, atlas_data, roi_ids)
% Purpose:
% Extract complete-ROI modulated gray matter volume in cubic millimeters
% from a CAT12 image on the fixed project atlas grid.
%
% Author: Qirui Zhang
% Created: 31 August 2026
%
%EXTRACT_REGIONAL_GMV Extract regional modulated gray matter volume.
%
% The atlas is already on the fixed GMV grid shared by both datasets. The
% complete retained ROI is used because structural voxel exclusion is not
% part of the approved workflow.

arguments
    gmv_path (1, 1) string
    atlas_data
    roi_ids (:, 1) double
end

gmv_info = niftiinfo(gmv_path);
gmv_data = double(niftiread(gmv_info));
gmv_data = squeeze(gmv_data);

if ~isequal(size(gmv_data), size(atlas_data))
    error('Block1:GMVSizeMismatch', ...
        'GMV image size does not match the atlas: %s', gmv_path);
end
voxel_volume_mm3 = prod(double(gmv_info.PixelDimensions(1:3)));
n_rois = numel(roi_ids);
gmv_mm3 = zeros(n_rois, 1);
total_roi_voxel_counts = zeros(n_rois, 1);
used_roi_voxel_counts = zeros(n_rois, 1);
excluded_roi_voxel_counts = zeros(n_rois, 1);

for roi_index = 1:n_rois
    total_roi_mask = atlas_data == roi_ids(roi_index);
    total_roi_voxel_counts(roi_index) = sum(total_roi_mask(:));
    used_roi_voxel_counts(roi_index) = total_roi_voxel_counts(roi_index);
    roi_values = gmv_data(total_roi_mask);
    gmv_mm3(roi_index) = sum(roi_values) * voxel_volume_mm3;
end
end

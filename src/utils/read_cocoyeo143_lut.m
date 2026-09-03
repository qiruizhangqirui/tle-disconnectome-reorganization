function atlas_table = read_cocoyeo143_lut(lut_path)
% Purpose:
% Read the cocoyeo143 lookup table and assign the encoded network label to
% each cortical or subcortical region.
%
% Author: Qirui Zhang
% Created: 31 August 2026
%
%READ_COCOYEO143_LUT Read and annotate the cocoyeo143 lookup table.
%
% lut_path identifies the approved cocoyeo143 text lookup table.
% atlas_table contains ROI_ID, ROI_Name, and Network in atlas order. Entries
% with nonpositive ROI IDs are removed. Cortical network labels are parsed
% from ROI names, and retained noncortical entries are labeled Subcortical.

arguments
    lut_path (1, 1) string
end

% Read and filter LUT entries.
lut_table = readtable(lut_path, ...
    'FileType', 'text', ...
    'Delimiter', ' ', ...
    'MultipleDelimsAsOne', true, ...
    'CommentStyle', '#', ...
    'ReadVariableNames', false);
roi_ids = lut_table.Var1;
roi_names = string(lut_table.Var2);
keep_rows = roi_ids > 0;
roi_ids = roi_ids(keep_rows);
roi_names = roi_names(keep_rows);

% Derive network labels from cortical ROI names.
network = repmat("Subcortical", size(roi_names));

is_cortical = startsWith(roi_names, "lh.") | startsWith(roi_names, "rh.");

cortical_names = roi_names(is_cortical);
network(is_cortical) = extractBefore( ...
    extractAfter(cortical_names, '7Networks_'), '_');

% Assemble the ordered atlas definition table.
atlas_table = table( ...
    roi_ids, roi_names, network, ...
    'VariableNames', {'ROI_ID', 'ROI_Name', 'Network'});
end

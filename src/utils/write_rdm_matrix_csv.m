function write_rdm_matrix_csv(output_path, matrix, roi_ids)
% Purpose:
% Write one square RDM with its ordered ROI IDs in the first row and column.
%
% Author: Qirui Zhang
% Created: 1 September 2026
%
%WRITE_RDM_MATRIX_CSV Write a square RDM with explicit ROI row and column IDs.

roi_ids = roi_ids(:);
output_directory = fileparts(output_path);
if ~isfolder(output_directory)
    mkdir(output_directory);
end
contents = cell(numel(roi_ids) + 1, numel(roi_ids) + 1);
contents(1, :) = cellstr(["ROI_ID", string(roi_ids')]);
contents(2:end, 1) = num2cell(roi_ids);
contents(2:end, 2:end) = num2cell(matrix);
writecell(contents, output_path);
end

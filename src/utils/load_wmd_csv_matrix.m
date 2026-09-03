function matrix = load_wmd_csv_matrix(csv_path)
% Purpose:
% Read one fixed cocoyeo143 NeMo connectivity CSV as a numeric matrix.
%
% Author: Qirui Zhang
% Created: 1 September 2026
%
%LOAD_WMD_CSV_MATRIX Load one full cocoyeo143 NeMo connectivity CSV.

contents = readcell(csv_path);
matrix = cell2mat(contents(2:end, 2:end));
end

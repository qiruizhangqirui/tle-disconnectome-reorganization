function plot_patient_rdm_montage(output_path, matrices, rdm_names, subject_id)
% Purpose:
% Export all patient RDM types in one diagnostic montage.
%
% Author: Qirui Zhang
% Created: 1 September 2026
%
%PLOT_PATIENT_RDM_MONTAGE Export all RDM types for one patient.

rdm_names = rdm_names(:);
output_directory = fileparts(output_path);
if ~isfolder(output_directory)
    mkdir(output_directory);
end

if numel(rdm_names) <= 4
    n_columns = numel(rdm_names);
elseif numel(rdm_names) <= 6
    n_columns = 3;
else
    n_columns = 4;
end
n_rows = ceil(numel(rdm_names) / n_columns);
figure_handle = figure('Visible', 'off', 'Color', 'white', ...
    'Position', [100, 100, 480 * n_columns, 450 * n_rows]);
layout = tiledlayout(figure_handle, n_rows, n_columns, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for rdm_index = 1:numel(rdm_names)
    matrix = matrices{rdm_index};
    axes_handle = nexttile(layout);
    imagesc(axes_handle, matrix);
    axis(axes_handle, 'image');
    axis(axes_handle, 'off');
    set(axes_handle, 'YDir', 'reverse');
    colormap(axes_handle, turbo(256));
    if any(matrix < 0, 'all')
        maximum_absolute_value = max(abs(matrix), [], 'all');
        clim(axes_handle, [-maximum_absolute_value, maximum_absolute_value]);
    elseif max(matrix, [], 'all') > min(matrix, [], 'all')
        clim(axes_handle, [min(matrix, [], 'all'), max(matrix, [], 'all')]);
    end
    colorbar(axes_handle);
    title(axes_handle, rdm_names(rdm_index), ...
        'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'normal');
end
title(layout, subject_id + " RDM Montage", ...
    'Interpreter', 'none', 'FontWeight', 'normal');
exportgraphics(figure_handle, output_path, 'Resolution', 250);
close(figure_handle);
end

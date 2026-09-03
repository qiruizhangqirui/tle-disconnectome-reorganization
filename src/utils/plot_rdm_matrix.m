function plot_rdm_matrix(output_path, matrix, roi_ids, plot_title)
% Purpose:
% Export one patient RDM as an individual matrix heatmap.
%
% Author: Qirui Zhang
% Created: 1 September 2026
%
%PLOT_RDM_MATRIX Export one high resolution RDM matrix figure.

roi_ids = roi_ids(:);
output_directory = fileparts(output_path);
if ~isfolder(output_directory)
    mkdir(output_directory);
end

figure_handle = figure('Visible', 'off', 'Color', 'white', ...
    'Position', [100, 100, 1400, 1200]);
axes_handle = axes(figure_handle);
imagesc(axes_handle, matrix);
axis(axes_handle, 'image');
set(axes_handle, 'YDir', 'reverse', 'FontSize', 9, 'TickLength', [0, 0]);
colormap(axes_handle, turbo(256));
if any(matrix < 0, 'all')
    maximum_absolute_value = max(abs(matrix), [], 'all');
    clim(axes_handle, [-maximum_absolute_value, maximum_absolute_value]);
elseif max(matrix, [], 'all') > min(matrix, [], 'all')
    clim(axes_handle, [min(matrix, [], 'all'), max(matrix, [], 'all')]);
end
n_retained = numel(roi_ids);
tick_positions = unique([1, round(linspace(20, n_retained, 5)), n_retained]);
xticks(axes_handle, tick_positions);
yticks(axes_handle, tick_positions);
xticklabels(axes_handle, string(tick_positions));
yticklabels(axes_handle, string(tick_positions));
xtickangle(axes_handle, 45);
xlabel(axes_handle, 'Retained ROIs');
ylabel(axes_handle, 'Retained ROIs');
title(axes_handle, plot_title, 'Interpreter', 'none', 'FontWeight', 'normal');
colorbar(axes_handle);
exportgraphics(figure_handle, output_path, 'Resolution', 300);
close(figure_handle);
end

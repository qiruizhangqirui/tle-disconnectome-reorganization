function plot_wmd_clinical_heatmap(effect_matrix, fdr_p_matrix, ...
    clinical_labels, analysis_labels, output_path)
% plot_wmd_clinical_heatmap: Summarize TJU clinical associations
%
% Purpose:
% Display the fixed 7 by 3 signed effect matrix with FDR-adjusted P values.
% Continuous effects are Spearman rho and categorical effects are signed r
% values derived from Welch t statistics.
%
% Author: Qirui Zhang
% Creation Date: 2026-08-31

% Preserve the analysis matrix order and construct a fixed diverging map.
row_positions = 1:7;
display_effects = effect_matrix;

blue = [0.230, 0.299, 0.754];
white = [1.000, 1.000, 1.000];
red = [0.706, 0.016, 0.150];
color_steps = 128;
diverging_map = [ ...
    linspace(blue(1), white(1), color_steps)', ...
    linspace(blue(2), white(2), color_steps)', ...
    linspace(blue(3), white(3), color_steps)'; ...
    linspace(white(1), red(1), color_steps)', ...
    linspace(white(2), red(2), color_steps)', ...
    linspace(white(3), red(3), color_steps)'];

% Draw the signed effects on a common scale from minus one to one.
figure('Color', 'w', 'Units', 'inches', 'Visible', 'off', ...
    'Position', [1, 1, 7.4, 5.2]);
axes_handle = axes('Position', [0.35, 0.16, 0.46, 0.72]);
image_handle = imagesc(axes_handle, display_effects);
image_handle.AlphaData = isfinite(display_effects);
set(axes_handle, 'Color', 'w');
colormap(axes_handle, diverging_map);
clim(axes_handle, [-1, 1]);
hold(axes_handle, 'on');

% Label every cell and outline associations with adjusted P below 0.05.
for clinical_index = 1:7
    display_row = row_positions(clinical_index);
    for analysis_index = 1:3
        effect_value = effect_matrix(clinical_index, analysis_index);
        fdr_p_value = fdr_p_matrix(clinical_index, analysis_index);
        fdr_p_text = string(compose('%.3g', fdr_p_value));
        if fdr_p_value < 0.05
            font_weight = 'bold';
        else
            font_weight = 'normal';
        end
        if abs(effect_value) > 0.55
            text_color = 'w';
        else
            text_color = 'k';
        end
        text(axes_handle, analysis_index, display_row, ...
            string(compose('%.2f', effect_value)) + newline + ...
            "P = " + fdr_p_text, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Color', text_color, 'FontName', 'Arial', ...
            'FontSize', 7.5, 'FontWeight', font_weight);
        if fdr_p_value < 0.05
            rectangle(axes_handle, 'Position', ...
                [analysis_index - 0.48, display_row - 0.48, 0.96, 0.96], ...
                'EdgeColor', 'k', 'LineWidth', 1.8);
        end
    end
end

% Apply the formal typography and export the 600 dpi figure.
xticks(axes_handle, 1:3);
xticklabels(axes_handle, analysis_labels);
xtickangle(axes_handle, 18);
yticks(axes_handle, row_positions);
yticklabels(axes_handle, clinical_labels);
set(axes_handle, 'FontName', 'Arial', 'FontSize', 8.5, ...
    'TickDir', 'out', 'Box', 'off');
axes_handle.Toolbar.Visible = 'off';
colorbar_handle = colorbar(axes_handle);
colorbar_handle.Label.String = 'Signed association effect';
title_obj = title(axes_handle, 'TJU clinical associations with WMD beta', 'FontSize', 10);
title_obj.Units = 'normalized';
title_obj.Position(2) = 1.05;
exportgraphics(gcf, output_path, 'Resolution', 600);
close(gcf);
end

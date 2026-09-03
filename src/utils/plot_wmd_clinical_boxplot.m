function plot_wmd_clinical_boxplot(beta_values, group_values, ...
    group_order, plot_title, plot_color, hedges_g, fdr_p_value, output_path)
% plot_wmd_clinical_boxplot: Plot a categorical clinical association
%
% Purpose:
% Plot WMD beta distributions in the fixed display order with individual
% observations, Hedges' g, and the FDR-adjusted P value.
%
% Author: Qirui Zhang
% Creation Date: 2026-08-31

% Standardize input orientation and create the fixed figure canvas.
beta_values = double(beta_values(:));
group_values = string(group_values(:));
group_order = string(group_order(:));

light_color = 0.55 + 0.45 * plot_color;
group_colors = [light_color; plot_color];
figure('Color', 'w', 'Units', 'inches', 'Visible', 'off', ...
    'Position', [1, 1, 7.4, 5.2]);
hold on;

% Draw each group in the fixed display order with deterministic jitter.
group_counts = zeros(2, 1);
for group_index = 1:2
    group_mask = group_values == group_order(group_index);
    group_counts(group_index) = sum(group_mask);
    group_beta = beta_values(group_mask);
    boxchart(group_index * ones(group_counts(group_index), 1), ...
        group_beta, 'BoxFaceColor', group_colors(group_index, :), ...
        'BoxFaceAlpha', 0.25, 'MarkerStyle', 'none', ...
        'LineWidth', 1.1);
    jitter = 0.16 * (2 * rand(group_counts(group_index), 1) - 1);
    scatter(group_index + jitter, group_beta, 30, ...
        group_colors(group_index, :), 'o', 'filled', ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
end
yline(0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1);

% Annotate the effect estimate, adjusted P value, and group sample sizes.
text(0.04, 0.96, "Hedges g = " + string(compose('%.3f', hedges_g)), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Arial', 'FontSize', 9);
fdr_p_text = string(compose('%.3g', fdr_p_value));
if fdr_p_value < 0.05
    p_font_weight = 'bold';
else
    p_font_weight = 'normal';
end
text(0.04, 0.89, "P = " + fdr_p_text, ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Arial', 'FontSize', 9, ...
    'FontWeight', p_font_weight);
text(0.04, 0.82, group_order(1) + " N = " + group_counts(1) + ...
    "; " + group_order(2) + " N = " + group_counts(2), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Arial', 'FontSize', 9);

% Apply the formal typography and export the 600 dpi figure.
xlim([0.5, 2.5]);
xticks([1, 2]);
xticklabels(cellstr(group_order));
ylabel('Standardized WMD beta');
title_obj = title(plot_title);
title_obj.Units = 'normalized';
title_obj.Position(2) = 1.04;
axes_handle = gca;
set(axes_handle, 'FontName', 'Arial', 'FontSize', 10, ...
    'LineWidth', 1, 'TickDir', 'out', 'Box', 'off');
axes_handle.Toolbar.Visible = 'off';
exportgraphics(gcf, output_path, 'Resolution', 600);
close(gcf);
end

function plot_wmd_clinical_scatter(x_values, beta_values, x_label, ...
    plot_title, plot_color, spearman_rho, fdr_p_value, output_path)
% plot_wmd_clinical_scatter: Plot a continuous clinical association
%
% Purpose:
% Plot complete clinical and WMD beta observations with a least-squares line,
% its 95% confidence band, Spearman rho, and the FDR-adjusted P value.
%
% Author: Qirui Zhang
% Creation Date: 2026-08-31

% Standardize input orientation and create the fixed figure canvas.
x_values = double(x_values(:));
beta_values = double(beta_values(:));

figure('Color', 'w', 'Units', 'inches', 'Visible', 'off', ...
    'Position', [1, 1, 7.4, 5.2]);
hold on;

% Draw the least-squares line, confidence band, observations, and zero line.
linear_model = fitlm(x_values, beta_values);
x_grid = linspace(min(x_values), max(x_values), 200)';
[predicted_beta, prediction_interval] = predict(linear_model, x_grid);
patch([x_grid; flipud(x_grid)], ...
    [prediction_interval(:, 1); flipud(prediction_interval(:, 2))], ...
    plot_color, 'FaceAlpha', 0.14, 'EdgeColor', 'none');
plot(x_grid, predicted_beta, 'Color', plot_color, 'LineWidth', 1.5);
scatter(x_values, beta_values, 34, plot_color, 'o', 'filled', ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
yline(0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1);

% Pad both axes around the complete plotted content.
x_range = max(x_values) - min(x_values);
x_padding = 0.06 * x_range;
xlim([min(x_values) - x_padding, max(x_values) + x_padding]);
y_content = [beta_values; prediction_interval(:); 0];
y_range = max(y_content) - min(y_content);
y_padding = 0.06 * y_range;
ylim([min(y_content) - y_padding, max(y_content) + y_padding]);

% Annotate the complete sample size, effect estimate, and adjusted P value.
text(0.04, 0.96, "N = " + string(numel(x_values)), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Arial', 'FontSize', 9);
text(0.04, 0.89, "Spearman rho = " + ...
    string(compose('%.3f', spearman_rho)), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Arial', 'FontSize', 9);
fdr_p_text = string(compose('%.3g', fdr_p_value));
if fdr_p_value < 0.05
    p_font_weight = 'bold';
else
    p_font_weight = 'normal';
end
text(0.04, 0.82, "P = " + fdr_p_text, ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Arial', 'FontSize', 9, ...
    'FontWeight', p_font_weight);

% Apply the formal typography and export the 600 dpi figure.
xlabel(x_label);
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

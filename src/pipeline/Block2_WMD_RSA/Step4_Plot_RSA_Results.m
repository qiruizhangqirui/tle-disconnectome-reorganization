%% Step4_Plot_RSA_Results: Create Publication Figures and Group Statistics for WMD RSA
%
% Purpose:
% Summarize center-specific whole-brain and pooled network WMD RSA results,
% then export group statistics and publication figures.
%
% Author: Qirui Zhang
% Created: 2026-08-31

%% Project paths and figure settings
% Define project-relative output paths and deterministic plotting settings.

block_directory = string(fileparts(mfilename('fullpath')));
project_root = string(fileparts(fileparts(fileparts(block_directory))));
result_root = fullfile(project_root, 'results', 'Block2_WMD_RSA');
figure_root = fullfile(result_root, 'figures');
if ~isfolder(figure_root)
    mkdir(figure_root);
end
distribution_directory = figure_root;
rng(12);

tju_color = [0.0000, 0.4470, 0.7410];
jlh_color = [0.8500, 0.3250, 0.0980];
point_size = 34;
jitter_width = 0.18;

%% Read formal RSA results
% Read whole-brain patient-level results for all analyses and both cohorts.

tju_post = readtable(fullfile(result_root, 'TJU', ...
    'POST_FC_RSA_Results.csv'), 'TextType', 'string');
tju_func_reorg = readtable(fullfile(result_root, 'TJU', ...
    'FUNC_REORG_RSA_Results.csv'), 'TextType', 'string');
tju_morph_reorg = readtable(fullfile(result_root, 'TJU', ...
    'MORPH_REORG_RSA_Results.csv'), 'TextType', 'string');
jlh_post = readtable(fullfile(result_root, 'JLH', ...
    'POST_FC_RSA_Results.csv'), 'TextType', 'string');
jlh_func_reorg = readtable(fullfile(result_root, 'JLH', ...
    'FUNC_REORG_RSA_Results.csv'), 'TextType', 'string');
jlh_morph_reorg = readtable(fullfile(result_root, 'JLH', ...
    'MORPH_REORG_RSA_Results.csv'), 'TextType', 'string');

%% Calculate one sample t tests
% Test each center-specific mean standardized WMD beta against zero.

analysis_names = ["POST_FC", "FUNC_REORG", "MORPH_REORG"];
whole_brain_cohort_names = ["TJU", "JLH"];
wmd_beta_vectors = { ...
    tju_post.Standardized_Beta_WMD, ...
    jlh_post.Standardized_Beta_WMD; ...
    tju_func_reorg.Standardized_Beta_WMD, ...
    jlh_func_reorg.Standardized_Beta_WMD; ...
    tju_morph_reorg.Standardized_Beta_WMD, ...
    jlh_morph_reorg.Standardized_Beta_WMD};

n_tests = numel(analysis_names) * numel(whole_brain_cohort_names);
test_analysis = strings(n_tests, 1);
test_cohort = strings(n_tests, 1);
test_n = zeros(n_tests, 1);
test_mean = zeros(n_tests, 1);
test_sd = zeros(n_tests, 1);
test_t = zeros(n_tests, 1);
test_df = zeros(n_tests, 1);
test_p = zeros(n_tests, 1);
test_d = zeros(n_tests, 1);
test_index = 0;
for analysis_index = 1:numel(analysis_names)
    for cohort_index = 1:numel(whole_brain_cohort_names)
        test_index = test_index + 1;
        values = wmd_beta_vectors{analysis_index, cohort_index};
        [~, ~, ~, statistics] = ttest(values, 0, 'Tail', 'both');
        p_value = two_sided_t_p_value( ...
            statistics.tstat, statistics.df);
        test_analysis(test_index) = analysis_names(analysis_index);
        test_cohort(test_index) = whole_brain_cohort_names(cohort_index);
        test_n(test_index) = numel(values);
        test_mean(test_index) = mean(values);
        test_sd(test_index) = std(values, 0);
        test_t(test_index) = statistics.tstat;
        test_df(test_index) = statistics.df;
        test_p(test_index) = p_value;
        test_d(test_index) = mean(values) / std(values, 0);
    end
end

test_table = table(test_analysis, test_cohort, test_n, test_mean, ...
    test_sd, test_t, test_df, test_p, test_d, ...
    'VariableNames', {'Analysis', 'Cohort', 'N', 'Mean_Beta_WMD', ...
    'SD_Beta_WMD', 'T_Statistic', 'Degrees_of_Freedom', ...
    'P_TwoSided', 'Cohens_D'});
writetable(test_table, fullfile(result_root, ...
    'Step4_OneSample_TTests.csv'));

%% Plot whole brain beta distributions by center
% Plot patient coefficients and center-specific two-sided P values.

whole_brain_tju_tables = {tju_post, tju_func_reorg, tju_morph_reorg};
whole_brain_jlh_tables = {jlh_post, jlh_func_reorg, jlh_morph_reorg};
whole_brain_titles = ["Postoperative FC RSA", ...
    "Functional reorganization dissimilarity RSA", ...
    "Morphometric reorganization dissimilarity RSA"];
cohort_colors = [tju_color; jlh_color];
for analysis_index = 1:numel(analysis_names)
    center_tables = {whole_brain_tju_tables{analysis_index}, ...
        whole_brain_jlh_tables{analysis_index}};
    figure('Color', 'w', 'Units', 'inches', ...
        'Position', [1, 1, 7.4, 5.2]);
    hold on;
    center_values = [];
    for center_index = 1:2
        center_table = center_tables{center_index};
        significant = center_table.Permutation_P_WMD < 0.05;
        x_values = center_index + jitter_width * ...
            (2 * rand(height(center_table), 1) - 1);
        boxchart(center_index * ones(height(center_table), 1), ...
            center_table.Standardized_Beta_WMD, ...
            'BoxFaceColor', cohort_colors(center_index, :), ...
            'BoxFaceAlpha', 0.18, 'MarkerStyle', 'none', ...
            'LineWidth', 1.1);
        scatter(x_values(~significant), ...
            center_table.Standardized_Beta_WMD(~significant), ...
            point_size, cohort_colors(center_index, :), ...
            'o', 'LineWidth', 0.9);
        scatter(x_values(significant), ...
            center_table.Standardized_Beta_WMD(significant), ...
            point_size, cohort_colors(center_index, :), 'o', 'filled', ...
            'MarkerEdgeColor', cohort_colors(center_index, :));
        center_values = [center_values; ...
            center_table.Standardized_Beta_WMD]; %#ok<AGROW>
    end
    yline(0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1);
    y_range = max(center_values) - min(center_values);
    ylim([min(center_values) - 0.10 * y_range, ...
        max(center_values) + 0.30 * y_range]);
    for center_index = 1:2
        p_value = test_table.P_TwoSided( ...
            test_table.Analysis == analysis_names(analysis_index) & ...
            test_table.Cohort == whole_brain_cohort_names(center_index));
        if p_value == 0
            p_label = "< 1e-15";
        elseif p_value < 0.01
            p_label = string(compose('%.2e', p_value));
        else
            p_label = string(compose('%.3g', p_value));
        end
        p_label = "P value" + newline + p_label;
        if p_value < 0.05
            p_font_weight = 'bold';
        else
            p_font_weight = 'normal';
        end
        text(center_index, max(center_values) + 0.13 * y_range, ...
            p_label, ...
            'HorizontalAlignment', 'center', 'FontName', 'Arial', ...
            'FontSize', 12, 'FontWeight', p_font_weight);
    end
    xlim([0.5, 2.5]);
    xticks([1, 2]);
    xticklabels({'TJU', 'JLH'});
    ylabel('Standardized WMD beta');
    title(whole_brain_titles(analysis_index));
    set(gca, 'FontName', 'Arial', 'FontSize', 10, ...
        'LineWidth', 1, 'TickDir', 'out', 'Box', 'off');
    exportgraphics(gcf, fullfile(distribution_directory, ...
        analysis_names(analysis_index) + ...
        "_Beta_By_Center.png"), 'Resolution', 600);
    close(gcf);

end

%% Plot postoperative FC unique delta R squared
% Plot PRE FC and WMD unique Delta R^2 separately within each cohort.

% Plot TJU distributions.
pre_r2_tju = tju_post.Unique_Delta_R2_PRE;
pre_p_sig_tju = tju_post.Permutation_P_PRE < 0.05;
[pre_r2_sorted_tju, pre_sort_order_tju] = sort(pre_r2_tju);
pre_p_sig_sorted_tju = pre_p_sig_tju(pre_sort_order_tju);

figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 7.4, 5.2]);
hold on;
scatter(find(~pre_p_sig_sorted_tju), pre_r2_sorted_tju(~pre_p_sig_sorted_tju), ...
    point_size, tju_color, 'o', 'LineWidth', 0.9);
scatter(find(pre_p_sig_sorted_tju), pre_r2_sorted_tju(pre_p_sig_sorted_tju), ...
    point_size, tju_color, 'o', 'filled', ...
    'MarkerEdgeColor', tju_color);
yline(0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1);
xlabel('Patients (sorted by PRE unique delta R^2)');
ylabel('PRE unique delta R^2');
set(gca, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1, ...
    'TickDir', 'out', 'Box', 'off');
exportgraphics(gcf, fullfile(figure_root, ...
    'POST_FC_PRE_Unique_Delta_R2_TJU.png'), 'Resolution', 600);
close(gcf);

wmd_r2_tju = tju_post.Unique_Delta_R2_WMD;
wmd_p_sig_tju = tju_post.Permutation_P_WMD < 0.05;
[wmd_r2_sorted_tju, wmd_sort_order_tju] = sort(wmd_r2_tju);
wmd_p_sig_sorted_tju = wmd_p_sig_tju(wmd_sort_order_tju);

figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 7.4, 5.2]);
hold on;
scatter(find(~wmd_p_sig_sorted_tju), wmd_r2_sorted_tju(~wmd_p_sig_sorted_tju), ...
    point_size, tju_color, 'o', 'LineWidth', 0.9);
scatter(find(wmd_p_sig_sorted_tju), wmd_r2_sorted_tju(wmd_p_sig_sorted_tju), ...
    point_size, tju_color, 'o', 'filled', ...
    'MarkerEdgeColor', tju_color);
yline(0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1);
xlabel('Patients (sorted by WMD unique delta R^2)');
ylabel('WMD unique delta R^2');
set(gca, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1, ...
    'TickDir', 'out', 'Box', 'off');
exportgraphics(gcf, fullfile(figure_root, ...
    'POST_FC_WMD_Unique_Delta_R2_TJU.png'), 'Resolution', 600);
close(gcf);

% Plot JLH distributions.
pre_r2_jlh = jlh_post.Unique_Delta_R2_PRE;
pre_p_sig_jlh = jlh_post.Permutation_P_PRE < 0.05;
[pre_r2_sorted_jlh, pre_sort_order_jlh] = sort(pre_r2_jlh);
pre_p_sig_sorted_jlh = pre_p_sig_jlh(pre_sort_order_jlh);

figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 7.4, 5.2]);
hold on;
scatter(find(~pre_p_sig_sorted_jlh), pre_r2_sorted_jlh(~pre_p_sig_sorted_jlh), ...
    point_size, jlh_color, 'o', 'LineWidth', 0.9);
scatter(find(pre_p_sig_sorted_jlh), pre_r2_sorted_jlh(pre_p_sig_sorted_jlh), ...
    point_size, jlh_color, 'o', 'filled', ...
    'MarkerEdgeColor', jlh_color);
yline(0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1);
xlabel('Patients (sorted by PRE unique delta R^2)');
ylabel('PRE unique delta R^2');
set(gca, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1, ...
    'TickDir', 'out', 'Box', 'off');
exportgraphics(gcf, fullfile(figure_root, ...
    'POST_FC_PRE_Unique_Delta_R2_JLH.png'), 'Resolution', 600);
close(gcf);

wmd_r2_jlh = jlh_post.Unique_Delta_R2_WMD;
wmd_p_sig_jlh = jlh_post.Permutation_P_WMD < 0.05;
[wmd_r2_sorted_jlh, wmd_sort_order_jlh] = sort(wmd_r2_jlh);
wmd_p_sig_sorted_jlh = wmd_p_sig_jlh(wmd_sort_order_jlh);

figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 7.4, 5.2]);
hold on;
scatter(find(~wmd_p_sig_sorted_jlh), wmd_r2_sorted_jlh(~wmd_p_sig_sorted_jlh), ...
    point_size, jlh_color, 'o', 'LineWidth', 0.9);
scatter(find(wmd_p_sig_sorted_jlh), wmd_r2_sorted_jlh(wmd_p_sig_sorted_jlh), ...
    point_size, jlh_color, 'o', 'filled', ...
    'MarkerEdgeColor', jlh_color);
yline(0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1);
xlabel('Patients (sorted by WMD unique delta R^2)');
ylabel('WMD unique delta R^2');
set(gca, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1, ...
    'TickDir', 'out', 'Box', 'off');
exportgraphics(gcf, fullfile(figure_root, ...
    'POST_FC_WMD_Unique_Delta_R2_JLH.png'), 'Resolution', 600);
close(gcf);


%% Read patient by network RSA results
% Read patient-level coefficients for the eight network-anchored models.

network_names = ["Visual", "Somatomotor", "DorsalAttention", ...
    "SalienceVentralAttention", "Limbic", "Control", ...
    "Default", "Subcortical"];
network_labels = {'Visual', 'Somatomotor', 'Dorsal attention', ...
    'Salience ventral attention', 'Limbic', 'Control', ...
    'Default', 'Subcortical'};
network_tju_tables = {
    readtable(fullfile(result_root, 'TJU', ...
        'POST_FC_Network_RSA_Results.csv'), 'TextType', 'string'), ...
    readtable(fullfile(result_root, 'TJU', ...
        'FUNC_REORG_Network_RSA_Results.csv'), 'TextType', 'string'), ...
    readtable(fullfile(result_root, 'TJU', ...
        'MORPH_REORG_Network_RSA_Results.csv'), 'TextType', 'string')};
network_jlh_tables = {
    readtable(fullfile(result_root, 'JLH', ...
        'POST_FC_Network_RSA_Results.csv'), 'TextType', 'string'), ...
    readtable(fullfile(result_root, 'JLH', ...
        'FUNC_REORG_Network_RSA_Results.csv'), 'TextType', 'string'), ...
    readtable(fullfile(result_root, 'JLH', ...
        'MORPH_REORG_Network_RSA_Results.csv'), 'TextType', 'string')};

%% Summarize network beta distributions
% Test each network mean and adjust P values within each eight-network family.

cohort_names = ["TJU", "JLH", "Combined"];
n_network_tests = numel(analysis_names) * numel(cohort_names) * ...
    numel(network_names);
summary_analysis = strings(n_network_tests, 1);
summary_cohort = strings(n_network_tests, 1);
summary_network = strings(n_network_tests, 1);
summary_network_order = zeros(n_network_tests, 1);
summary_n = zeros(n_network_tests, 1);
summary_mean = zeros(n_network_tests, 1);
summary_absolute_mean = zeros(n_network_tests, 1);
summary_sd = zeros(n_network_tests, 1);
summary_se = zeros(n_network_tests, 1);
summary_ci_lower = zeros(n_network_tests, 1);
summary_ci_upper = zeros(n_network_tests, 1);
summary_t = zeros(n_network_tests, 1);
summary_df = zeros(n_network_tests, 1);
summary_p = zeros(n_network_tests, 1);
summary_fdr_p = zeros(n_network_tests, 1);
summary_d = zeros(n_network_tests, 1);
summary_rank = zeros(n_network_tests, 1);
summary_direction = strings(n_network_tests, 1);
summary_index = 0;

for analysis_index = 1:numel(analysis_names)
    for cohort_index = 1:numel(cohort_names)
        if cohort_names(cohort_index) == "TJU"
            cohort_table = network_tju_tables{analysis_index};
        elseif cohort_names(cohort_index) == "JLH"
            cohort_table = network_jlh_tables{analysis_index};
        else
            cohort_table = [network_tju_tables{analysis_index}; ...
                network_jlh_tables{analysis_index}];
        end
        family_rows = zeros(numel(network_names), 1);
        for network_index = 1:numel(network_names)
            summary_index = summary_index + 1;
            family_rows(network_index) = summary_index;
            values = cohort_table.Standardized_Beta_WMD( ...
                cohort_table.Network == network_names(network_index));
            [~, ~, confidence_interval, statistics] = ttest(values);
            p_value = two_sided_t_p_value( ...
                statistics.tstat, statistics.df);
            summary_analysis(summary_index) = analysis_names(analysis_index);
            summary_cohort(summary_index) = cohort_names(cohort_index);
            summary_network(summary_index) = network_names(network_index);
            summary_network_order(summary_index) = network_index;
            summary_n(summary_index) = numel(values);
            summary_mean(summary_index) = mean(values);
            summary_absolute_mean(summary_index) = abs(mean(values));
            summary_sd(summary_index) = std(values, 0);
            summary_se(summary_index) = std(values, 0) / sqrt(numel(values));
            summary_ci_lower(summary_index) = confidence_interval(1);
            summary_ci_upper(summary_index) = confidence_interval(2);
            summary_t(summary_index) = statistics.tstat;
            summary_df(summary_index) = statistics.df;
            summary_p(summary_index) = p_value;
            summary_d(summary_index) = mean(values) / std(values, 0);
            if mean(values) > 0
                summary_direction(summary_index) = "Positive";
            elseif mean(values) < 0
                summary_direction(summary_index) = "Negative";
            else
                summary_direction(summary_index) = "Zero";
            end
        end

        family_p = summary_p(family_rows);
        summary_fdr_p(family_rows) = mafdr( ...
            family_p, 'BHFDR', true);

        [~, importance_order] = sort( ...
            summary_absolute_mean(family_rows), 'descend');
        family_rank = zeros(numel(network_names), 1);
        family_rank(importance_order) = 1:numel(network_names);
        summary_rank(family_rows) = family_rank;
    end
end

network_summary = table(summary_analysis, summary_cohort, ...
    summary_network, summary_network_order, summary_n, summary_mean, ...
    summary_absolute_mean, summary_sd, summary_se, summary_ci_lower, ...
    summary_ci_upper, summary_t, summary_df, summary_p, summary_fdr_p, ...
    summary_d, summary_rank, summary_direction, ...
    'VariableNames', {'Analysis', 'Cohort', 'Network', 'Network_Order', ...
    'N', 'Mean_Beta_WMD', 'Absolute_Mean_Beta_WMD', 'SD_Beta_WMD', ...
    'SE_Beta_WMD', 'CI95_Lower', 'CI95_Upper', 'T_Statistic', ...
    'Degrees_of_Freedom', 'P_TwoSided', ...
    'FDR_P_Across_8_Networks', 'Cohens_D', 'Importance_Rank', ...
    'Direction'});
writetable(network_summary, fullfile(result_root, ...
    'Network_Group_Summary.csv'));

%% Plot network beta distributions
% Plot center-specific patient distributions with FDR-adjusted P values.

analysis_titles = ["Postoperative FC RSA", ...
    "FC change profile RSA", "GMV change RSA"];
combined_color = [0.25, 0.25, 0.25];
for analysis_index = 1:numel(analysis_names)
    center_tables = {network_tju_tables{analysis_index}, ...
        network_jlh_tables{analysis_index}};
    for center_index = 1:2
            center_table = center_tables{center_index};
            figure('Color', 'w', 'Units', 'inches', ...
                'Position', [1, 1, 7.4, 5.2]);
            hold on;
            boxchart(center_table.Network_Order, ...
                center_table.Standardized_Beta_WMD, ...
                'BoxFaceColor', cohort_colors(center_index, :), ...
                'BoxFaceAlpha', 0.18, 'MarkerStyle', 'none', ...
                'LineWidth', 1.0);
            x_values = center_table.Network_Order + 0.13 * ...
                (2 * rand(height(center_table), 1) - 1);
            scatter(x_values, center_table.Standardized_Beta_WMD, ...
                18, cohort_colors(center_index, :), 'o', 'filled');
            yline(0, ':', 'Color', [0.35, 0.35, 0.35]);
            xticks(1:8);
            xticklabels(network_labels);
            xtickangle(35);
            ylabel('Standardized WMD beta');
            title(analysis_titles(analysis_index) + " " + ...
                cohort_names(center_index));
            set(gca, 'FontName', 'Arial', 'FontSize', 9, ...
                'LineWidth', 1, 'TickDir', 'out', 'Box', 'off');
            center_summary = network_summary( ...
                network_summary.Analysis == analysis_names(analysis_index) & ...
                network_summary.Cohort == cohort_names(center_index), :);
            y_limits = ylim;
            y_range = range(y_limits);
            annotation_y = y_limits(2) + 0.11 * y_range;
            ylim([y_limits(1), y_limits(2) + 0.30 * y_range]);
            for network_index = 1:8
                adjusted_p = center_summary. ...
                    FDR_P_Across_8_Networks(network_index);
                if adjusted_p == 0
                    p_label = "< 1e-15";
                elseif adjusted_p < 0.01
                    p_label = string(compose('%.2e', adjusted_p));
                else
                    p_label = string(compose('%.3g', adjusted_p));
                end
                p_label = "P value" + newline + p_label;
                if adjusted_p < 0.05
                    p_font_weight = 'bold';
                else
                    p_font_weight = 'normal';
                end
                text(network_index, annotation_y, p_label, ...
                    'HorizontalAlignment', 'center', ...
                    'Color', cohort_colors(center_index, :), ...
                    'FontName', 'Arial', 'FontSize', 10, ...
                    'FontWeight', p_font_weight);
            end
            exportgraphics(gcf, fullfile(figure_root, ...
                analysis_names(analysis_index) + "_Network_Beta_" + ...
                cohort_names(center_index) + ".png"), 'Resolution', 600);
            close(gcf);
    end

end

%% Plot network forest plots
% Plot network means, 95% CIs, and FDR-adjusted P values by cohort.

forest_colors = [0.122, 0.467, 0.706; 0.851, 0.325, 0.098; ...
    0.25, 0.25, 0.25];
for analysis_index = 1:numel(analysis_names)
    for cohort_index = 1:3
        figure('Color', 'w', 'Units', 'inches', ...
            'Position', [1, 1, 7.4, 5.2]);
        hold on;
        rows = network_summary( ...
            network_summary.Analysis == analysis_names(analysis_index) & ...
            network_summary.Cohort == cohort_names(cohort_index), :);
        y_values = (1:8)';
        lower_error = rows.Mean_Beta_WMD - rows.CI95_Lower;
        upper_error = rows.CI95_Upper - rows.Mean_Beta_WMD;
        errorbar(rows.Mean_Beta_WMD, ...
            y_values, lower_error, upper_error, 'horizontal', ...
            'LineStyle', 'none', 'Marker', 'o', ...
            'Color', forest_colors(cohort_index, :), ...
            'MarkerFaceColor', forest_colors(cohort_index, :), ...
            'LineWidth', 1.0, 'MarkerSize', 4);
        xline(0, ':', 'Color', [0.35, 0.35, 0.35]);
        data_x_limits = xlim;
        data_x_ticks = xticks;
        data_x_range = range(data_x_limits);
        p_column_x = data_x_limits(2) + 0.18 * data_x_range;
        xlim([data_x_limits(1), ...
            data_x_limits(2) + 0.36 * data_x_range]);
        xticks(data_x_ticks);
        ylim([0.00, 8.70]);
        yticks(1:8);
        yticklabels(network_labels);
        text(p_column_x, 0.52, 'P value', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontName', 'Arial', ...
            'FontSize', 11, 'FontWeight', 'bold');
        for network_index = 1:8
            adjusted_p = rows.FDR_P_Across_8_Networks(network_index);
            if adjusted_p == 0
                p_label = "< 1e-15";
            elseif adjusted_p < 0.01
                p_label = string(compose('%.2e', adjusted_p));
            else
                p_label = string(compose('%.3g', adjusted_p));
            end
            if adjusted_p < 0.05
                p_font_weight = 'bold';
            else
                p_font_weight = 'normal';
            end
            text(p_column_x, network_index, p_label, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontName', 'Arial', ...
                'FontSize', 10, 'FontWeight', p_font_weight);
        end
        set(gca, 'YDir', 'reverse', 'FontName', 'Arial', 'FontSize', 11, ...
            'LineWidth', 1, 'TickDir', 'out', 'Box', 'off');
        xlabel('Mean standardized WMD beta and 95% CI', 'FontSize', 12);
        title(analysis_titles(analysis_index) + " " + ...
            cohort_names(cohort_index), 'FontSize', 12.5, 'FontWeight', 'bold');
        exportgraphics(gcf, fullfile(figure_root, ...
            analysis_names(analysis_index) + "_Network_Forest_" + ...
            cohort_names(cohort_index) + ".png"), 'Resolution', 600);
        close(gcf);
    end
end

function p_value = two_sided_t_p_value(t_statistic, degrees_of_freedom)
% Compute a numerically stable two-sided Student t P value.
p_value = betainc(degrees_of_freedom / ...
    (degrees_of_freedom + t_statistic^2), ...
    degrees_of_freedom / 2, 0.5);
end

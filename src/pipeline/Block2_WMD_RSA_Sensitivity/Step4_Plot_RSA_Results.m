%% Step4_Plot_RSA_Results: Create Publication Figures and Group Statistics for Spatial Sensitivity RSA
%
% Purpose:
% Summarize center-specific whole-brain RSA results for the Euclidean distance
% and cavity distance difference sensitivity branches.
%
% Author: Qirui Zhang
% Created: 2026-08-31

%% Configuration and paths
% Define isolated branch outputs and deterministic plotting settings.

block_directory = string(fileparts(mfilename('fullpath')));
project_root = string(fileparts(fileparts(fileparts(block_directory))));
result_root = fullfile(project_root, 'results', 'Block2_WMD_RSA_Sensitivity');

branches = ["Euclidean_Distance", "Cavity_Distance_Difference"];

tju_color = [0.0000, 0.4470, 0.7410];
jlh_color = [0.8500, 0.3250, 0.0980];
cohort_colors = [tju_color; jlh_color];
point_size = 34;
jitter_width = 0.18;

analysis_names = ["POST_FC", "FUNC_REORG", "MORPH_REORG"];
analysis_titles = [ ...
    "Postoperative FC Sensitivity RSA", ...
    "Functional Reorganization Dissimilarity Sensitivity RSA", ...
    "Morphometric Reorganization Dissimilarity Sensitivity RSA"];
cohort_names = ["TJU", "JLH"];

%% Run Analysis for Each Sensitivity Branch
% Process each sensitivity covariate in a separate result branch.

for branch_index = 1:numel(branches)
    branch_name = branches(branch_index);
    branch_dir = fullfile(result_root, branch_name);
    figure_root = fullfile(branch_dir, 'figures');
    if ~isfolder(figure_root)
        mkdir(figure_root);
    end

    fprintf('\n========================================\n');
    fprintf('Processing Branch: %s\n', branch_name);
    fprintf('========================================\n');

    % Read patient-level results for all analyses and both cohorts.
    tju_post = readtable(fullfile(branch_dir, 'TJU', 'POST_FC_RSA_Results.csv'), 'TextType', 'string');
    tju_func = readtable(fullfile(branch_dir, 'TJU', 'FUNC_REORG_RSA_Results.csv'), 'TextType', 'string');
    tju_morph = readtable(fullfile(branch_dir, 'TJU', 'MORPH_REORG_RSA_Results.csv'), 'TextType', 'string');

    jlh_post = readtable(fullfile(branch_dir, 'JLH', 'POST_FC_RSA_Results.csv'), 'TextType', 'string');
    jlh_func = readtable(fullfile(branch_dir, 'JLH', 'FUNC_REORG_RSA_Results.csv'), 'TextType', 'string');
    jlh_morph = readtable(fullfile(branch_dir, 'JLH', 'MORPH_REORG_RSA_Results.csv'), 'TextType', 'string');

    tables_tju = {tju_post, tju_func, tju_morph};
    tables_jlh = {jlh_post, jlh_func, jlh_morph};

    %% 1. Group-Level One-Sample t-Tests for WMD and Sensitivity Covariate
    % Test center-specific mean standardized coefficients against zero.
    test_rows = {};
    for a_idx = 1:numel(analysis_names)
        a_name = analysis_names(a_idx);
        tab_tju = tables_tju{a_idx};
        tab_jlh = tables_jlh{a_idx};

        if a_name == "POST_FC"
            pred_list = ["PRE", "WMD", "Sensitivity_Covariate"];
            var_tju = {tab_tju.Standardized_Beta_PRE, tab_tju.Standardized_Beta_WMD, tab_tju.Standardized_Beta_Sensitivity_Covariate};
            var_jlh = {tab_jlh.Standardized_Beta_PRE, tab_jlh.Standardized_Beta_WMD, tab_jlh.Standardized_Beta_Sensitivity_Covariate};
        else
            pred_list = ["WMD", "Sensitivity_Covariate"];
            var_tju = {tab_tju.Standardized_Beta_WMD, tab_tju.Standardized_Beta_Sensitivity_Covariate};
            var_jlh = {tab_jlh.Standardized_Beta_WMD, tab_jlh.Standardized_Beta_Sensitivity_Covariate};
        end

        for p_idx = 1:numel(pred_list)
            p_name = pred_list(p_idx);
            v_tju = var_tju{p_idx};
            v_jlh = var_jlh{p_idx};

            cohort_vecs = {v_tju, v_jlh};
            for c_idx = 1:numel(cohort_names)
                c_name = cohort_names(c_idx);
                vals = cohort_vecs{c_idx};
                [~, ~, ~, stats] = ttest(vals, 0, 'Tail', 'both');
                p_val = two_sided_t_p_value(stats.tstat, stats.df);
                test_rows = [test_rows; {a_name, p_name, c_name, numel(vals), ...
                    mean(vals), std(vals, 0), stats.tstat, stats.df, p_val, mean(vals) / std(vals, 0)}]; %#ok<AGROW>
            end
        end
    end

    test_table = cell2table(test_rows, 'VariableNames', ...
        {'Analysis', 'Predictor', 'Cohort', 'N', 'Mean_Beta', 'SD_Beta', ...
         'T_Statistic', 'Degrees_of_Freedom', 'P_TwoSided', 'Cohens_D'});
    ttest_path = fullfile(figure_root, 'Step4_OneSample_TTests.csv');
    writetable(test_table, ttest_path);
    fprintf('Exported: %s\n', ttest_path);

    %% 2. Whole-Brain Beta Distributions by Center for WMD and Covariate
    % Plot patient coefficients and center-specific two-sided P values.
    rng(12);
    beta_plot_specs = struct();
    beta_plot_specs(1).pred_name = "WMD";
    beta_plot_specs(1).beta_var = "Standardized_Beta_WMD";
    beta_plot_specs(1).p_var = "Permutation_P_WMD";
    beta_plot_specs(1).file_tag = "";
    beta_plot_specs(1).ylabel_str = "Standardized WMD beta";

    beta_plot_specs(2).pred_name = "Sensitivity_Covariate";
    beta_plot_specs(2).beta_var = "Standardized_Beta_Sensitivity_Covariate";
    beta_plot_specs(2).p_var = "Permutation_P_Sensitivity_Covariate";
    beta_plot_specs(2).file_tag = "_COV";
    beta_plot_specs(2).ylabel_str = "Standardized Covariate beta";

    for a_idx = 1:numel(analysis_names)
        a_name = analysis_names(a_idx);
        tab_tju = tables_tju{a_idx};
        tab_jlh = tables_jlh{a_idx};
        center_tabs = {tab_tju, tab_jlh};

        for bp_idx = 1:numel(beta_plot_specs)
            spec = beta_plot_specs(bp_idx);
            pred_n = spec.pred_name;
            b_var = spec.beta_var;
            p_var = spec.p_var;
            f_tag = spec.file_tag;
            y_lbl = spec.ylabel_str;

            % Plot TJU and JLH distributions in one panel.
            figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 7.4, 5.2]);
            hold on;
            center_vals = [];
            for c_idx = 1:2
                tab = center_tabs{c_idx};
                sig = tab.(p_var) < 0.05;
                vals_c = tab.(b_var);
                x_pts = c_idx + jitter_width * (2 * rand(height(tab), 1) - 1);
                boxchart(c_idx * ones(height(tab), 1), vals_c, ...
                    'BoxFaceColor', cohort_colors(c_idx, :), 'BoxFaceAlpha', 0.18, ...
                    'MarkerStyle', 'none', 'LineWidth', 1.1);
                scatter(x_pts(~sig), vals_c(~sig), ...
                    point_size, cohort_colors(c_idx, :), 'o', 'LineWidth', 0.9);
                scatter(x_pts(sig), vals_c(sig), ...
                    point_size, cohort_colors(c_idx, :), 'o', 'filled', ...
                    'MarkerEdgeColor', cohort_colors(c_idx, :));
                center_vals = [center_vals; vals_c]; %#ok<AGROW>
            end
            yline(0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1);
            y_rng = max(center_vals) - min(center_vals);
            ylim([min(center_vals) - 0.10 * y_rng, max(center_vals) + 0.30 * y_rng]);

            for c_idx = 1:2
                sub_t = test_table(test_table.Analysis == a_name & ...
                    test_table.Predictor == pred_n & ...
                    test_table.Cohort == cohort_names(c_idx), :);
                p_val = sub_t.P_TwoSided(1);
                if p_val == 0
                    p_lbl = "< 1e-15";
                elseif p_val < 0.01
                    p_lbl = string(compose('%.2e', p_val));
                else
                    p_lbl = string(compose('%.3g', p_val));
                end
                p_lbl = "P value" + newline + p_lbl;
                if p_val < 0.05
                    p_weight = 'bold';
                else
                    p_weight = 'normal';
                end
                text(c_idx, max(center_vals) + 0.13 * y_rng, p_lbl, ...
                    'HorizontalAlignment', 'center', 'FontName', 'Arial', ...
                    'FontSize', 12, 'FontWeight', p_weight);
            end
            xlim([0.5, 2.5]);
            xticks([1, 2]);
            xticklabels({'TJU', 'JLH'});
            ylabel(y_lbl);
            set(gca, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1, ...
                'TickDir', 'out', 'Box', 'off', 'Position', [0.14, 0.13, 0.78, 0.73]);
            
            if pred_n == "WMD"
                t_str = analysis_titles(a_idx) + " (" + strrep(branch_name, '_', ' ') + ")";
            else
                t_str = analysis_titles(a_idx) + " - Covariate (" + strrep(branch_name, '_', ' ') + ")";
            end
            th = title(t_str);
            set(th, 'Units', 'normalized', 'Position', [0.5, 1.07, 0]);
            exportgraphics(gcf, fullfile(figure_root, a_name + f_tag + "_Beta_By_Center.png"), 'Resolution', 600);
            close(gcf);
        end
    end

    %% 3. Variance Partitioning (Combined Unique Delta R^2 Plots Across Centers)
    % Plot pooled patient-level unique Delta R^2 with cohort-specific colors.
    var_items = struct();
    
    % Postoperative FC model: PRE FC unique variance.
    var_items(1).fname = "POST_FC_PRE_Unique_Delta_R2";
    var_items(1).r2_tju = tju_post.Unique_Delta_R2_PRE;
    var_items(1).p_tju = tju_post.Permutation_P_PRE;
    var_items(1).r2_jlh = jlh_post.Unique_Delta_R2_PRE;
    var_items(1).p_jlh = jlh_post.Permutation_P_PRE;
    var_items(1).lbl = "PRE unique delta R^2";
    var_items(1).model_name = "Postoperative FC RSA";

    % Postoperative FC model: WMD unique variance.
    var_items(2).fname = "POST_FC_WMD_Unique_Delta_R2";
    var_items(2).r2_tju = tju_post.Unique_Delta_R2_WMD;
    var_items(2).p_tju = tju_post.Permutation_P_WMD;
    var_items(2).r2_jlh = jlh_post.Unique_Delta_R2_WMD;
    var_items(2).p_jlh = jlh_post.Permutation_P_WMD;
    var_items(2).lbl = "WMD unique delta R^2";
    var_items(2).model_name = "Postoperative FC RSA";

    % Postoperative FC model: sensitivity covariate unique variance.
    var_items(3).fname = "POST_FC_COV_Unique_Delta_R2";
    var_items(3).r2_tju = tju_post.Unique_Delta_R2_Sensitivity_Covariate;
    var_items(3).p_tju = tju_post.Permutation_P_Sensitivity_Covariate;
    var_items(3).r2_jlh = jlh_post.Unique_Delta_R2_Sensitivity_Covariate;
    var_items(3).p_jlh = jlh_post.Permutation_P_Sensitivity_Covariate;
    var_items(3).lbl = "Covariate unique delta R^2";
    var_items(3).model_name = "Postoperative FC RSA";

    % Functional reorganization model: WMD unique variance.
    var_items(4).fname = "FUNC_REORG_WMD_Unique_Delta_R2";
    var_items(4).r2_tju = tju_func.Unique_Delta_R2_WMD;
    var_items(4).p_tju = tju_func.Permutation_P_WMD;
    var_items(4).r2_jlh = jlh_func.Unique_Delta_R2_WMD;
    var_items(4).p_jlh = jlh_func.Permutation_P_WMD;
    var_items(4).lbl = "WMD unique delta R^2";
    var_items(4).model_name = "Functional Reorganization RSA";

    % Functional reorganization model: sensitivity covariate unique variance.
    var_items(5).fname = "FUNC_REORG_COV_Unique_Delta_R2";
    var_items(5).r2_tju = tju_func.Unique_Delta_R2_Sensitivity_Covariate;
    var_items(5).p_tju = tju_func.Permutation_P_Sensitivity_Covariate;
    var_items(5).r2_jlh = jlh_func.Unique_Delta_R2_Sensitivity_Covariate;
    var_items(5).p_jlh = jlh_func.Permutation_P_Sensitivity_Covariate;
    var_items(5).lbl = "Covariate unique delta R^2";
    var_items(5).model_name = "Functional Reorganization RSA";

    % Morphometric reorganization model: WMD unique variance.
    var_items(6).fname = "MORPH_REORG_WMD_Unique_Delta_R2";
    var_items(6).r2_tju = tju_morph.Unique_Delta_R2_WMD;
    var_items(6).p_tju = tju_morph.Permutation_P_WMD;
    var_items(6).r2_jlh = jlh_morph.Unique_Delta_R2_WMD;
    var_items(6).p_jlh = jlh_morph.Permutation_P_WMD;
    var_items(6).lbl = "WMD unique delta R^2";
    var_items(6).model_name = "Morphometric Reorganization RSA";

    % Morphometric reorganization model: sensitivity covariate unique variance.
    var_items(7).fname = "MORPH_REORG_COV_Unique_Delta_R2";
    var_items(7).r2_tju = tju_morph.Unique_Delta_R2_Sensitivity_Covariate;
    var_items(7).p_tju = tju_morph.Permutation_P_Sensitivity_Covariate;
    var_items(7).r2_jlh = jlh_morph.Unique_Delta_R2_Sensitivity_Covariate;
    var_items(7).p_jlh = jlh_morph.Permutation_P_Sensitivity_Covariate;
    var_items(7).lbl = "Covariate unique delta R^2";
    var_items(7).model_name = "Morphometric Reorganization RSA";

    for v_idx = 1:numel(var_items)
        fname = var_items(v_idx).fname;
        r2_tju = var_items(v_idx).r2_tju;
        p_tju = var_items(v_idx).p_tju;
        r2_jlh = var_items(v_idx).r2_jlh;
        p_jlh = var_items(v_idx).p_jlh;
        lbl = var_items(v_idx).lbl;
        model_name = var_items(v_idx).model_name;

        r2_all = [r2_tju; r2_jlh];
        p_all = [p_tju; p_jlh];
        cohort_id = [repmat("TJU", numel(r2_tju), 1); repmat("JLH", numel(r2_jlh), 1)];

        [r2_sorted, sort_order] = sort(r2_all);
        p_sorted = p_all(sort_order);
        cohort_sorted = cohort_id(sort_order);

        is_tju = cohort_sorted == "TJU";
        is_jlh = cohort_sorted == "JLH";
        is_sig = p_sorted < 0.05;

        figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 7.8, 5.2]);
        hold on;

        idx_tju_nonsig = find(is_tju & ~is_sig);
        idx_tju_sig = find(is_tju & is_sig);
        scatter(idx_tju_nonsig, r2_sorted(idx_tju_nonsig), ...
            point_size, tju_color, 'o', 'LineWidth', 0.9);
        scatter(idx_tju_sig, r2_sorted(idx_tju_sig), ...
            point_size, tju_color, 'o', 'filled', 'MarkerEdgeColor', tju_color);

        idx_jlh_nonsig = find(is_jlh & ~is_sig);
        idx_jlh_sig = find(is_jlh & is_sig);
        scatter(idx_jlh_nonsig, r2_sorted(idx_jlh_nonsig), ...
            point_size, jlh_color, 'o', 'LineWidth', 0.9);
        scatter(idx_jlh_sig, r2_sorted(idx_jlh_sig), ...
            point_size, jlh_color, 'o', 'filled', 'MarkerEdgeColor', jlh_color);

        yline(0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1);
        xlabel(sprintf('Patients (N=88, sorted by %s)', lbl));
        ylabel(lbl);

        dummy_tju = scatter(nan, nan, point_size, tju_color, 'o', 'filled', 'MarkerEdgeColor', tju_color);
        dummy_jlh = scatter(nan, nan, point_size, jlh_color, 'o', 'filled', 'MarkerEdgeColor', jlh_color);
        legend([dummy_tju, dummy_jlh], {'TJU (N=52)', 'JLH (N=36)'}, ...
            'Location', 'northwest', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 9);

        set(gca, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1, ...
            'TickDir', 'out', 'Box', 'off', 'Position', [0.14, 0.13, 0.78, 0.73]);

        plot_title_str = sprintf('%s: %s (%s)', model_name, lbl, strrep(branch_name, '_', ' '));
        th = title(plot_title_str, 'FontName', 'Arial', 'FontSize', 10.5, 'FontWeight', 'normal');
        set(th, 'Units', 'normalized', 'Position', [0.5, 1.07, 0]);

        exportgraphics(gcf, fullfile(figure_root, fname + ".png"), 'Resolution', 600);
        close(gcf);
    end
end

fprintf('\nAll sensitivity plotting and group statistics completed successfully!\n');

function p_value = two_sided_t_p_value(t_statistic, degrees_of_freedom)
% Compute a numerically stable two-sided Student t P value.
p_value = betainc(degrees_of_freedom / ...
    (degrees_of_freedom + t_statistic^2), ...
    degrees_of_freedom / 2, 0.5);
end

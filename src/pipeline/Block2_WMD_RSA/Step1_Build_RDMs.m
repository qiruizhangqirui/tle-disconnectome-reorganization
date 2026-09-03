%% Step1_Build_RDMs: Build Patient RDMs for Primary WMD RSA
%
% Purpose:
% Construct five patient-specific Representational Dissimilarity Matrices
% (RDMs) for the TJU and JLH cohorts using the same scientific sequence.
% The analysis excludes the 27 cerebellar regions and each patient's
% surgery-affected regions with at least 5 percent cavity overlap.
%
% 1. PRE_FC_RDM: One minus preoperative FC.
% 2. POST_FC_RDM: One minus postoperative FC.
% 3. FUNC_REORG_RDM: One minus the Pearson correlation between regional FC
%    change profiles after excluding the two regions defining each pair.
% 4. MORPH_REORG_RDM: Absolute pairwise difference in signed longitudinal
%    nodal morphometric reorganization percentage.
% 5. WMD_RDM: Continuous NeMo ChaCo disconnection among retained regions.
%
% Author: Qirui Zhang
% Created: 1 September 2026

%% Configuration

export_rdm_figures = true; % Export one image per RDM and one montage per patient.

%% Project paths and fixed cohort configuration

block_directory = string(fileparts(mfilename('fullpath')));
project_root = string(fileparts(fileparts(fileparts(block_directory))));
addpath(fullfile(project_root, 'src', 'utils'));

dataset_names = ["TJU"; "JLH"];
metadata_sheet_names = ["TJU_metadata"; "JLH_metadata"];
nemo_filename_suffixes = [ ...
    "_nemo_output_ifod2act_chacoconn_cocoyeo143subj_mean.csv"; ...
    "_POST1_Run1_ROI_nemo_output_ifod2act_chacoconn_cocoyeo143subj_mean.csv"];

%% Atlas definition

atlas_table = read_cocoyeo143_lut( ...
    fullfile(project_root, 'Data', 'atlas', 'cocoyeo143_LUT.txt'));
full_atlas_roi_ids = atlas_table.ROI_ID; % Ordered 143-region atlas axis.
noncerebellar_mask = ~startsWith(atlas_table.ROI_Name, "Cbm");
analysis_roi_ids = full_atlas_roi_ids(noncerebellar_mask); % Ordered 116-region analysis axis.

%% Dataset and patient RDM construction

for dataset_index = 1:numel(dataset_names)
    dataset = dataset_names(dataset_index);
    fprintf('\nDataset %s\n', dataset);

    block1_root = fullfile(project_root, 'results', ...
        'Block1_TLE_surgery_reorganization', dataset, 'Original');
    output_root = fullfile(project_root, 'results', 'Block2_WMD_RSA', dataset);
    rdm_output_directory = fullfile(output_root, 'intermediate');
    individual_figure_root = fullfile(output_root, 'figures', 'individual');
    montage_figure_directory = fullfile(output_root, 'figures', 'patient_montage');

    metadata = readtable( ...
        fullfile(project_root, 'Data', 'metadata', dataset + "_metadata.xlsx"), ...
        'Sheet', metadata_sheet_names(dataset_index), ...
        'VariableNamingRule', 'preserve');
    subject_ids = strtrim(string(metadata.("Subject ID")));

    % The Block 1 overlap table is the fixed source of patient-specific ROI exclusion.
    surgery_overlap_table = readtable( ...
        fullfile(block1_root, 'Surgery_ROI_Overlap.csv'), ...
        'VariableNamingRule', 'preserve', 'TextType', 'string');

    for subject_index = 1:numel(subject_ids)
        subject_id = subject_ids(subject_index);
        fprintf('Patient %02d/%02d: %s\n', ...
            subject_index, numel(subject_ids), subject_id);

        %% Load and align the patient inputs

        % The Block 1 MAT file supplies FC, GMV, and their shared 116-region axis.
        block2_input = load( ...
            fullfile(block1_root, 'block2_inputs', ...
            subject_id + "_Block2_Input.mat"), ...
            'roi_ids', 'pre_raw_fc', 'post_raw_fc', ...
            'pre_gmv_mm3', 'post_gmv_mm3');
        subject_overlap = surgery_overlap_table( ...
            string(surgery_overlap_table.Subject_ID) == subject_id, :);
        [~, overlap_order] = ismember( ...
            block2_input.roi_ids, subject_overlap.ROI_ID);
        surgery_affected = logical( ...
            subject_overlap.Surgery_Affected(overlap_order));
        analysis_rows = ismember(block2_input.roi_ids, analysis_roi_ids);
        retained_rows = analysis_rows & ~surgery_affected;
        retained_roi_ids = block2_input.roi_ids(retained_rows);
        surgery_affected_roi_ids = ...
            block2_input.roi_ids(analysis_rows & surgery_affected);
        [~, retained_roi_indices] = ismember( ...
            retained_roi_ids, analysis_roi_ids); % Positions on the 116-region analysis axis.

        retained_pre_fc = ...
            block2_input.pre_raw_fc(retained_rows, retained_rows);
        retained_post_fc = ...
            block2_input.post_raw_fc(retained_rows, retained_rows);
        retained_pre_gmv = block2_input.pre_gmv_mm3(retained_rows);
        retained_post_gmv = block2_input.post_gmv_mm3(retained_rows);
        n_rois = numel(retained_roi_ids);

        %% Construct the five primary RDMs

        pre_fc_rdm = 1 - retained_pre_fc; % Pairwise preoperative FC dissimilarity.
        pre_fc_rdm(1:n_rois + 1:end) = 0;

        post_fc_rdm = 1 - retained_post_fc; % Pairwise postoperative FC dissimilarity.
        post_fc_rdm(1:n_rois + 1:end) = 0;

        % Compare regional longitudinal FC change profiles without the direct pair.
        delta_fc_profile = retained_post_fc - retained_pre_fc;
        func_reorg_rdm = zeros(n_rois, n_rois);
        for first_roi = 1:n_rois - 1
            for second_roi = first_roi + 1:n_rois
                profile_mask = true(n_rois, 1);
                profile_mask([first_roi, second_roi]) = false;
                profile_correlation = corr( ...
                    delta_fc_profile(first_roi, profile_mask)', ...
                    delta_fc_profile(second_roi, profile_mask)');
                func_reorg_rdm(first_roi, second_roi) = ...
                    1 - profile_correlation;
                func_reorg_rdm(second_roi, first_roi) = ...
                    1 - profile_correlation;
            end
        end

        % NMR is signed POST minus PRE GMV change expressed as a percentage of PRE.
        nmr_vector = 100 * (retained_post_gmv(:) - retained_pre_gmv(:)) ./ ...
            retained_pre_gmv(:);
        morph_reorg_rdm = abs(nmr_vector - nmr_vector');
        morph_reorg_rdm(1:n_rois + 1:end) = 0;

        % ChaCo values are used directly as continuous pairwise WMD dissimilarities.
        nemo_filename = subject_id + nemo_filename_suffixes(dataset_index);
        full_wmd = load_wmd_csv_matrix(fullfile( ...
            project_root, 'Data', dataset, 'nemo', nemo_filename));
        noncerebellar_wmd = full_wmd(noncerebellar_mask, noncerebellar_mask);
        wmd_rdm = noncerebellar_wmd( ...
            retained_roi_indices, retained_roi_indices);
        wmd_rdm(1:n_rois + 1:end) = 0;

        matrices = {pre_fc_rdm, post_fc_rdm, func_reorg_rdm, ...
            morph_reorg_rdm, wmd_rdm};
        matrix_names = ["PRE_FC_RDM", "POST_FC_RDM", ...
            "FUNC_REORG_RDM", "MORPH_REORG_RDM", "WMD_RDM"];

        %% Write patient RDMs, metadata, and approved diagnostic figures

        for matrix_index = 1:numel(matrix_names)
            matrix_filename = subject_id + "_" + ...
                matrix_names(matrix_index) + ".csv";
            write_rdm_matrix_csv( ...
                fullfile(rdm_output_directory, matrix_filename), ...
                matrices{matrix_index}, retained_roi_ids);
            if export_rdm_figures
                plot_rdm_matrix(fullfile(individual_figure_root, ...
                    matrix_names(matrix_index), ...
                    subject_id + "_" + matrix_names(matrix_index) + ".png"), ...
                    matrices{matrix_index}, retained_roi_ids, ...
                    subject_id + " " + matrix_names(matrix_index));
            end
        end

        metadata_row = table( ...
            dataset, subject_id, n_rois, ...
            strjoin(string(retained_roi_ids'), ';'), ...
            strjoin(string(surgery_affected_roi_ids(:)'), ';'), ...
            'VariableNames', {'Dataset', 'Subject_ID', 'N_ROIs', ...
            'Retained_ROI_IDs', 'Excluded_ROI_IDs'});
        for matrix_index = 1:numel(matrix_names)
            column_name = char(matrix_names(matrix_index) + "_File");
            metadata_row.(column_name) = subject_id + "_" + ...
                matrix_names(matrix_index) + ".csv";
        end
        writetable(metadata_row, fullfile(rdm_output_directory, ...
            subject_id + "_RDM_Metadata.csv"));

        if export_rdm_figures
            plot_patient_rdm_montage(fullfile(montage_figure_directory, ...
                subject_id + "_RDM_Montage.png"), ...
                matrices, matrix_names, subject_id);
        end
    end
end

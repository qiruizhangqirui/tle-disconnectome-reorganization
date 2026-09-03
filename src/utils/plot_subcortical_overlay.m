function [a, cb] = plot_subcortical_overlay(subcortical_values, subcortical_mask, varargin)
% plot_subcortical_overlay
%
% Purpose:
% Render a dual-layer subcortical mesh with a cavity mask and colorbar.
%
% Author: Qirui Zhang
% Created: 31 August 2026
%
% Usage:
%   [a, cb] = plot_subcortical_overlay(subcortical_values, subcortical_mask, varargin)
%
% Inputs:
%   subcortical_values - Vector of 14 values (excluding ventricles) or 16 values.
%   subcortical_mask   - Logical vector of 14 (or 16) elements indicating whether each
%                        subcortical structure is in the group-level cavity mask.
%
% Optional Name-Value Pairs:
%   ventricles            - 'False' [default] or 'True'.
%   color_range           - [min, max] range for statistical overlay.
%   cmap                  - Colormap name ('RdBu_r' [default], 'Reds', 'Blues').
%   mask_color            - RGB color for masked structures (default [0.28, 0.28, 0.28]).
%   mesh_background_color - RGB color for neutral background (default [0.94, 0.94, 0.94]).
%   background            - Figure background color (default 'white').
%   label_text            - Colorbar title string (default "").

p = inputParser;
addParameter(p, 'ventricles', 'False', @ischar);
addParameter(p, 'color_range', [], @isnumeric);
addParameter(p, 'cmap', 'RdBu_r', @ischar);
addParameter(p, 'mask_color', [0.28, 0.28, 0.28], @isnumeric);
addParameter(p, 'mesh_background_color', [0.94, 0.94, 0.94], @isnumeric);
addParameter(p, 'background', 'white', @ischar);
addParameter(p, 'label_text', "", @(x) ischar(x) || isstring(x));

parse(p, varargin{:});
in = p.Results;

subcortical_values = double(subcortical_values(:));
subcortical_mask = logical(subcortical_mask(:));

% Load subcortical meshes
surf_lh = SurfStatReadSurf('sctx_lh');
surf_rh = SurfStatReadSurf('sctx_rh');

% Expand to 16 structures (including ventricles)
counts = [867; 1419; 3012; 3784; 1446; 4003; 3726; 7653; ...
          838; 1457; 3208; 3742; 1373; 3871; 3699; 7180];

if strcmp(in.ventricles, 'False')
    if numel(subcortical_values) ~= 14 || numel(subcortical_mask) ~= 14
        error('plot_subcortical_overlay:InvalidInputSize', ...
            'Expected 14 values when ventricles is False.');
    end
    data16 = nan(16, 1);
    data16([1:7, 9:15]) = subcortical_values;
    mask16 = false(16, 1);
    mask16([1:7, 9:15]) = subcortical_mask;
else
    if numel(subcortical_values) ~= 16 || numel(subcortical_mask) ~= 16
        error('plot_subcortical_overlay:InvalidInputSize', ...
            'Expected 16 values when ventricles is True.');
    end
    data16 = subcortical_values;
    mask16 = subcortical_mask;
end

% Assign values and masks to vertices
vert_values_cell = cell(16, 1);
vert_mask_cell = cell(16, 1);
for i = 1:16
    vert_values_cell{i} = repmat(data16(i), counts(i), 1);
    vert_mask_cell{i} = repmat(mask16(i), counts(i), 1);
end
vert_values = vertcat(vert_values_cell{:});
vert_mask = vertcat(vert_mask_cell{:});

% Determine color range
if isempty(in.color_range)
    valid_unmasked = vert_values(~vert_mask & isfinite(vert_values));
    if isempty(valid_unmasked)
        color_range = [-1, 1];
    else
        color_range = [min(valid_unmasked), max(valid_unmasked)];
        if color_range(1) == color_range(2)
            color_range = color_range(1) + [-1, 1];
        end
    end
else
    color_range = double(in.color_range);
end

% Load the approved ENIGMA colormap
cmap_matrix = enigma_colormap(in.cmap);
n_cmap_colors = size(cmap_matrix, 1);

% Build RGB matrix for all vertices
n_verts = length(vert_values);
rgb_matrix = repmat(in.mesh_background_color, n_verts, 1);

% Layer 1: Dark gray mask for surgery-affected structures
rgb_matrix(vert_mask, :) = repmat(in.mask_color, sum(vert_mask), 1);

% Layer 2: Overlay statistical values for unmasked structures with finite values
stat_indices = find(~vert_mask & isfinite(vert_values));
if ~isempty(stat_indices)
    vals = vert_values(stat_indices);
    norm_vals = (vals - color_range(1)) / (color_range(2) - color_range(1));
    norm_vals = min(max(norm_vals, 0), 1);
    color_idx = round(norm_vals * (n_cmap_colors - 1)) + 1;
    color_idx = min(max(color_idx, 1), n_cmap_colors);
    rgb_matrix(stat_indices, :) = cmap_matrix(color_idx, :);
end

% When ventricles is False, skip ventricles by setting vertex RGB to NaN
if strcmp(in.ventricles, 'False')
    ventricle_indices = [18258:25910, 44099:length(vert_values)];
    rgb_matrix(ventricle_indices, :) = NaN;
end

vl = 1:size(surf_lh.coord, 2);
vr = (1:size(surf_rh.coord, 2)) + max(size(surf_lh.coord, 2));

% Tight, compact layout with narrow spacing between the 4 views
h = 0.58;
w = 0.235;
y_pos = 0.28;
x_pos = [0.02, 0.26, 0.51, 0.75];

% Axis 1: Left Lateral
a(1) = axes('Position', [x_pos(1), y_pos, w, h]);
trisurf(surf_lh.tri, surf_lh.coord(1, :), surf_lh.coord(2, :), surf_lh.coord(3, :), ...
    'FaceVertexCData', rgb_matrix(vl, :), 'FaceColor', 'interp', 'EdgeColor', 'none');
view(-90, 0);
daspect([1, 1, 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;

% Axis 2: Left Medial
a(2) = axes('Position', [x_pos(2), y_pos, w, h]);
trisurf(surf_lh.tri, surf_lh.coord(1, :), surf_lh.coord(2, :), surf_lh.coord(3, :), ...
    'FaceVertexCData', rgb_matrix(vl, :), 'FaceColor', 'interp', 'EdgeColor', 'none');
view(90, 0);
daspect([1, 1, 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;

% Axis 3: Right Lateral
a(3) = axes('Position', [x_pos(3), y_pos, w, h]);
trisurf(surf_rh.tri, surf_rh.coord(1, :), surf_rh.coord(2, :), surf_rh.coord(3, :), ...
    'FaceVertexCData', rgb_matrix(vr, :), 'FaceColor', 'interp', 'EdgeColor', 'none');
view(-90, 0);
daspect([1, 1, 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;

% Axis 4: Right Medial
a(4) = axes('Position', [x_pos(4), y_pos, w, h]);
trisurf(surf_rh.tri, surf_rh.coord(1, :), surf_rh.coord(2, :), surf_rh.coord(3, :), ...
    'FaceVertexCData', rgb_matrix(vr, :), 'FaceColor', 'interp', 'EdgeColor', 'none');
view(90, 0);
daspect([1, 1, 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;

for i = 1:4
    set(a(i), 'Tag', ['SurfStatView ' num2str(i)]);
    set(a(i), 'CLim', color_range);
    colormap(a(i), cmap_matrix);
end

% Add Colorbar
cb = colorbar(a(4), 'location', 'South');
set(cb, 'Position', [0.35, 0.10, 0.30, 0.04]);
set(cb, 'XAxisLocation', 'bottom');
set(cb, 'FontName', 'Arial', 'FontSize', 11);
h_title = get(cb, 'Title');
set(h_title, 'String', in.label_text, 'FontName', 'Arial', 'FontSize', 11);

set(gcf, 'Color', in.background, 'InvertHardcopy', 'off');

end

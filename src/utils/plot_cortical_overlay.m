function [a, cb] = plot_cortical_overlay(cortical_values, cortical_mask, varargin)
% plot_cortical_overlay
%
% Purpose:
% Render a dual-layer cortical surface with a cavity mask and colorbar.
%
% Author: Qirui Zhang
% Created: 31 August 2026
%
% Usage:
%   [a, cb] = plot_cortical_overlay(cortical_values, cortical_mask, varargin)
%
% Inputs:
%   cortical_values - Vector of 100 values in ENIGMA Schaefer-100 order.
%   cortical_mask   - Logical vector of 100 elements indicating whether each
%                     parcel is in the group-level cavity mask.
%
% Optional Name-Value Pairs:
%   surface_name          - Surface mesh ('conte69' [default], 'fsa5').
%   color_range           - [min, max] range for statistical overlay.
%   cmap                  - Colormap name ('RdBu_r' [default], 'Reds', 'Blues').
%   mask_color            - RGB color for masked parcels (default [0.28, 0.28, 0.28]).
%   mesh_background_color - RGB color for neutral background (default [0.94, 0.94, 0.94]).
%   background            - Figure background color (default 'white').
%   label_text            - Colorbar title string (default "").

p = inputParser;
addParameter(p, 'surface_name', 'conte69', @ischar);
addParameter(p, 'color_range', [], @isnumeric);
addParameter(p, 'cmap', 'RdBu_r', @ischar);
addParameter(p, 'mask_color', [0.28, 0.28, 0.28], @isnumeric);
addParameter(p, 'mesh_background_color', [0.94, 0.94, 0.94], @isnumeric);
addParameter(p, 'background', 'white', @ischar);
addParameter(p, 'label_text', "", @(x) ischar(x) || isstring(x));

parse(p, varargin{:});
in = p.Results;

cortical_values = double(cortical_values(:));
cortical_mask = logical(cortical_mask(:));

if numel(cortical_values) ~= 100 || numel(cortical_mask) ~= 100
    error('plot_cortical_overlay:InvalidInputSize', ...
        'Expected 100 Schaefer parcels for cortical values and mask.');
end

% Load surface mesh
if strcmp(in.surface_name, 'fsa5')
    surf = SurfStatAvSurf({'fsa5_lh', 'fsa5_rh'});
    mapping_name = 'schaefer_100_fsa5';
elseif strcmp(in.surface_name, 'conte69')
    surf = SurfStatAvSurf({'conte69_lh', 'conte69_rh'});
    mapping_name = 'schaefer_100_conte69';
else
    error('plot_cortical_overlay:UnsupportedSurface', ...
        'Unsupported surface name: %s.', in.surface_name);
end

% Map parcel values and mask to surface vertices
vert_values = parcel_to_surface(cortical_values, mapping_name, nan);
vert_mask = parcel_to_surface(double(cortical_mask), mapping_name, 0) >= 0.5;

% Determine color range
if isempty(in.color_range)
    valid_unmasked_vals = vert_values(~vert_mask & isfinite(vert_values));
    if isempty(valid_unmasked_vals)
        color_range = [-1, 1];
    else
        color_range = [min(valid_unmasked_vals), max(valid_unmasked_vals)];
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

% Layer 1: Apply dark charcoal gray mask for surgery-affected parcels
rgb_matrix(vert_mask, :) = repmat(in.mask_color, sum(vert_mask), 1);

% Layer 2: Overlay statistical values for unmasked parcels with finite values
stat_vert_indices = find(~vert_mask & isfinite(vert_values));
if ~isempty(stat_vert_indices)
    vals = vert_values(stat_vert_indices);
    norm_vals = (vals - color_range(1)) / (color_range(2) - color_range(1));
    norm_vals = min(max(norm_vals, 0), 1);
    color_idx = round(norm_vals * (n_cmap_colors - 1)) + 1;
    color_idx = min(max(color_idx, 1), n_cmap_colors);
    rgb_matrix(stat_vert_indices, :) = cmap_matrix(color_idx, :);
end

% Set up rendering parameters
vl = 1:(n_verts / 2);
vr = vl + (n_verts / 2);
t = size(surf.tri, 1);
tl = 1:(t / 2);
tr = tl + (t / 2);

% Tight, compact layout with narrow spacing between the 4 views
h = 0.58;
w = 0.235;
y_pos = 0.28;
x_pos = [0.02, 0.26, 0.51, 0.75];

% Axis 1: Left Lateral
a(1) = axes('Position', [x_pos(1), y_pos, w, h]);
trisurf(surf.tri(tl, :), surf.coord(1, vl), surf.coord(2, vl), surf.coord(3, vl), ...
    'FaceVertexCData', rgb_matrix(vl, :), 'FaceColor', 'interp', 'EdgeColor', 'none');
view(-90, 0);
daspect([1, 1, 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;

% Axis 2: Left Medial
a(2) = axes('Position', [x_pos(2), y_pos, w, h]);
trisurf(surf.tri(tl, :), surf.coord(1, vl), surf.coord(2, vl), surf.coord(3, vl), ...
    'FaceVertexCData', rgb_matrix(vl, :), 'FaceColor', 'interp', 'EdgeColor', 'none');
view(90, 0);
daspect([1, 1, 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;

% Axis 3: Right Lateral
a(3) = axes('Position', [x_pos(3), y_pos, w, h]);
trisurf(surf.tri(tr, :) - (n_verts / 2), surf.coord(1, vr), surf.coord(2, vr), surf.coord(3, vr), ...
    'FaceVertexCData', rgb_matrix(vr, :), 'FaceColor', 'interp', 'EdgeColor', 'none');
view(-90, 0);
daspect([1, 1, 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;

% Axis 4: Right Medial
a(4) = axes('Position', [x_pos(4), y_pos, w, h]);
trisurf(surf.tri(tr, :) - (n_verts / 2), surf.coord(1, vr), surf.coord(2, vr), surf.coord(3, vr), ...
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

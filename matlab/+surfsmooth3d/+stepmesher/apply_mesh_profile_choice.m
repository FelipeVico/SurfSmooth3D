function [opts, choice] = apply_mesh_profile_choice(opts, choiceName)
%APPLY_MESH_PROFILE_CHOICE Apply a documented scaffold meshing preset.
%
% Usage:
%   opts = surfsmooth3d.stepmesher.apply_mesh_profile_choice(opts, 'curvature_balanced');
%
% See also surfsmooth3d.stepmesher.mesh_profile_choices.

if nargin < 1 || isempty(opts)
    opts = struct();
end
if nargin < 2 || isempty(choiceName)
    choiceName = 'rigid';
end

choices = surfsmooth3d.stepmesher.mesh_profile_choices();
names = lower(string(choices.name));
match = find(names == lower(string(choiceName)), 1);
if isempty(match)
    error('stepmesher:meshProfileChoice:unknownChoice', ...
        'Unknown mesh profile choice "%s". Use surfsmooth3d.stepmesher.mesh_profile_choices() to list valid choices.', ...
        string(choiceName));
end

choice = choices(match, :);
opts.mesh_profile = char(choice.mesh_profile);
opts.curvature_points = choice.curvature_points;
opts.hmin_mode = char(choice.hmin_mode);
opts.hmin_fraction = choice.hmin_fraction;
opts.hmin_edge_scale = choice.hmin_edge_scale;
opts.cad_edge_tol_fraction = choice.cad_edge_tol_fraction;
opts.hmax_fraction = choice.hmax_fraction;
opts.mesh_size_extend_from_boundary = ...
    choice.mesh_size_extend_from_boundary;
opts.gmsh_algorithm = char(choice.gmsh_algorithm);
opts.optimize = choice.optimize;
opts.min_angle_deg = choice.min_angle_deg;
opts.max_edge_ratio = choice.max_edge_ratio;
opts.enforce_quality = choice.enforce_quality;
end

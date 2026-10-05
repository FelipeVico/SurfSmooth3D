function opts = default_options(opts)
%DEFAULT_OPTIONS Fill option defaults copied from the Python reference.

defaults = struct();
defaults.order = 8;
defaults.mesh_fraction = 0.10;
defaults.refinement_level = 0;
defaults.edge_gll_order = 16;
defaults.mesh_profile = 'rigid';
defaults.curvature_points = 50;
defaults.hmin_mode = 'fraction';
defaults.hmin_fraction = 0.01;
defaults.hmin_edge_scale = 0.25;
defaults.cad_edge_tol_fraction = 1e-6;
defaults.hmax_fraction = 0.10;
defaults.mesh_size_extend_from_boundary = true;
defaults.gmsh_algorithm = 'frontal_delaunay';
defaults.optimize = true;
defaults.min_angle_deg = 15;
defaults.max_edge_ratio = 6;
defaults.enforce_quality = false;
defaults.occt_precision = 1e-6;
defaults.occt_maxprecision = 1e-6;
privateDir = fileparts(mfilename('fullpath'));
occInfoExe = fullfile(privateDir, 'step_mesher_occ_info');
if isfile(occInfoExe)
    defaults.occ_info_exe = occInfoExe;
else
    defaults.occ_info_exe = '';
end
defaults.sameparameter = true;
defaults.surfacecurve_mode = '3d_preferred';
defaults.verbose = false;

names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(opts, name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end
end

%% STEP Mesher MATLAB Driver
% MATLAB analog of the original Python notebook entry point.
%
% This script assumes:
%   1. step_mesher_mex has been built by CMake.
%   2. SurfSmooth3D/matlab is available for geometry and basis functions.

clearvars
close all
clc

%% Locate repositories
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));

addpath(fullfile(repoRoot, 'matlab'));
setup_surfsmooth3d();

%% User parameters
% Choose one STEP file by uncommenting the desired line.
% These examples are copied from the Python step_mesher repo, plus
% torus_ellipse.step from the release test set.
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_geometry_1.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_geometry_2.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_geometry_manas.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'UV_curved_bezier.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_curved_uv.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'cylinder_bezier_v2.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'intersection_two_balls.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'piecewise_bezier_cap.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'spheres_intersect.step');
stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'torus_ellipse.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'Multi_res.step');



% Choose one scaffold/adaptivity launch preset.
% To try curvature-based Gmsh sizing, change this to 'curvature_balanced'.
% Other useful choices are 'curvature_coarse',
% 'curvature_balanced_local', 'curvature_balanced_local_cadedge',
% 'curvature_balanced_wide', 'curvature_fine', and 'curvature_guarded'.
% The full menu is printed when the driver starts.
MESH_PROFILE_CHOICE = 'rigid';
% MESH_PROFILE_CHOICE = 'curvature_balanced';
% MESH_PROFILE_CHOICE = 'curvature_balanced_local';
% MESH_PROFILE_CHOICE = 'curvature_balanced_local_cadedge';
% MESH_PROFILE_CHOICE = 'curvature_balanced_wide';
% MESH_PROFILE_CHOICE = 'curvature_fine';
% MESH_PROFILE_CHOICE = 'curvature_guarded';

% Optional curvature-profile overrides. Leave [] to use the selected
% preset defaults. For 'curvature_balanced_local', the defaults are:
%   curvature_points = 50
%   hmin_mode = 'fraction'
%   hmin_fraction = 0.01
%   hmin_edge_scale = 0.25
%   cad_edge_tol_fraction = 1e-6
%   hmax_fraction = 0.10
%   mesh_size_extend_from_boundary = false
CURVATURE_POINTS_OVERRIDE = [];
HMIN_MODE_OVERRIDE = '';
HMIN_FRACTION_OVERRIDE = [];
HMIN_EDGE_SCALE_OVERRIDE = [];
CAD_EDGE_TOL_FRACTION_OVERRIDE = [];
HMAX_FRACTION_OVERRIDE = [];
MESH_SIZE_EXTEND_FROM_BOUNDARY_OVERRIDE = [];
% CURVATURE_POINTS_OVERRIDE = 35;
% HMIN_MODE_OVERRIDE = 'cad_edge';
% HMIN_FRACTION_OVERRIDE = 0.01;
% HMIN_EDGE_SCALE_OVERRIDE = 0.25;
% CAD_EDGE_TOL_FRACTION_OVERRIDE = 1e-6;
% HMAX_FRACTION_OVERRIDE = 0.12;
% MESH_SIZE_EXTEND_FROM_BOUNDARY_OVERRIDE = false;

opts = struct();
opts.order = 8;
opts.mesh_fraction = 0.050;
opts.refinement_level = 0;
opts.edge_gll_order = opts.order;
[opts, meshProfileChoice] = surfsmooth3d.stepmesher.apply_mesh_profile_choice( ...
    opts, MESH_PROFILE_CHOICE);
if ~isempty(CURVATURE_POINTS_OVERRIDE)
    opts.curvature_points = CURVATURE_POINTS_OVERRIDE;
end
if ~isempty(HMIN_MODE_OVERRIDE)
    opts.hmin_mode = HMIN_MODE_OVERRIDE;
end
if ~isempty(HMIN_FRACTION_OVERRIDE)
    opts.hmin_fraction = HMIN_FRACTION_OVERRIDE;
end
if ~isempty(HMIN_EDGE_SCALE_OVERRIDE)
    opts.hmin_edge_scale = HMIN_EDGE_SCALE_OVERRIDE;
end
if ~isempty(CAD_EDGE_TOL_FRACTION_OVERRIDE)
    opts.cad_edge_tol_fraction = CAD_EDGE_TOL_FRACTION_OVERRIDE;
end
if ~isempty(HMAX_FRACTION_OVERRIDE)
    opts.hmax_fraction = HMAX_FRACTION_OVERRIDE;
end
if ~isempty(MESH_SIZE_EXTEND_FROM_BOUNDARY_OVERRIDE)
    opts.mesh_size_extend_from_boundary = ...
        MESH_SIZE_EXTEND_FROM_BOUNDARY_OVERRIDE;
end
opts.occt_precision = 1e-6;
opts.occt_maxprecision = 1e-6;
opts.sameparameter = true;
opts.surfacecurve_mode = '3d_preferred';
opts.verbose = false;

stepStem = erase_file_extension(get_filename(stepFile));
outDir = fullfile(repoRoot, 'outputs', 'step', stepStem);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
xyzwFile = fullfile(outDir, sprintf('%s_order%d_xyzw.txt', stepStem, opts.order));
packageName = sprintf('%s_order%d', stepStem, opts.order);
packageDir = fullfile(outDir, packageName);

%% Run mesher
fprintf('Available scaffold/adaptivity choices:\n');
disp(surfsmooth3d.stepmesher.mesh_profile_choices());
fprintf('Selected choice: %s\n', meshProfileChoice.name);
fprintf('  description: %s\n', meshProfileChoice.description);
fprintf('  C++ profile: %s\n', meshProfileChoice.mesh_profile);
if strcmp(meshProfileChoice.mesh_profile, "curvature")
    fprintf(['  curvature_points=%d, hmin_mode=%s, ', ...
        'hmin_fraction=%.4g, hmin_edge_scale=%.4g, ', ...
        'cad_edge_tol_fraction=%.4g, hmax_fraction=%.4g\n'], ...
        opts.curvature_points, opts.hmin_mode, opts.hmin_fraction, ...
        opts.hmin_edge_scale, opts.cad_edge_tol_fraction, ...
        opts.hmax_fraction);
    fprintf(['  extend_from_boundary=%d, gmsh_algorithm=%s, ', ...
        'optimize=%d\n'], opts.mesh_size_extend_from_boundary, ...
        opts.gmsh_algorithm, opts.optimize);
    fprintf(['  min_angle=%.4g, max_edge_ratio=%.4g, ', ...
        'enforce_quality=%d\n'], opts.min_angle_deg, ...
        opts.max_edge_ratio, opts.enforce_quality);
end
fprintf('\n');

fprintf('STEP file: %s\n', stepFile);
fprintf('Order:     %d\n', opts.order);
fprintf('Mesh frac: %.6g\n', opts.mesh_fraction);
fprintf('Profile:   %s (%s)\n\n', opts.mesh_profile, MESH_PROFILE_CHOICE);

mesh = surfsmooth3d.stepmesher.mesh_step(stepFile, opts);

fprintf('Mesh complete.\n');
fprintf('  patches:          %d\n', mesh.npatches);
fprintf('  nodes per patch:  %d\n', mesh.nodes_per_patch);
fprintf('  total nodes:      %d\n', size(mesh.xyz, 2));
fprintf('  CAD area:         %.16e\n', mesh.meta.cad_area);
fprintf('  hmin/hmax eff:    %.16e / %.16e\n', ...
    mesh.meta.hmin_effective, mesh.meta.hmax_effective);
fprintf('  CAD edge min:     raw %.16e, meaningful %.16e (%d/%d)\n', ...
    mesh.meta.min_raw_cad_edge_length, ...
    mesh.meta.min_meaningful_cad_edge_length, ...
    mesh.meta.meaningful_cad_edge_count, mesh.meta.cad_edge_count);
fprintf('  max proj dist:    %.16e\n', mesh.meta.max_projection_distance);
fprintf('  mean proj dist:   %.16e\n\n', mesh.meta.mean_projection_distance);
fprintf('  min tri angle:    %.6e deg\n', mesh.meta.min_scaffold_angle_deg);
fprintf('  max edge ratio:   %.6e\n', mesh.meta.max_scaffold_edge_ratio);
fprintf('  bad angle/ratio:  %d / %d\n\n', ...
    mesh.meta.bad_min_angle_count, mesh.meta.bad_edge_ratio_count);

%% Convert to fmm3dbie formats
[srcvals, norders, iptype, meta] = surfsmooth3d.stepmesher.to_srcvals(mesh); %#ok<ASGLU>
S = surfsmooth3d.stepmesher.to_surfer(mesh);

fprintf('surfer constructed.\n');
fprintf('  S.npatches:       %d\n', S.npatches);
fprintf('  S.npts:           %d\n\n', S.npts);

%% Validation and export
stats = surfsmooth3d.stepmesher.validate_area_flux(mesh);
payload_xyzw = surfsmooth3d.stepmesher.export_xyzw(mesh, xyzwFile);
payload_package = surfsmooth3d.stepmesher.export_geom_package(mesh, packageDir, packageName);

fprintf('Validation:\n');
fprintf('  quadrature area:  %.16e\n', stats.area);
fprintf('  CAD area:         %.16e\n', stats.cad_area);
fprintf('  rel area error:   %.16e\n', stats.relative_area_error);
fprintf('  normal flux:      [%.6e, %.6e, %.6e]\n\n', ...
    stats.normal_flux(1), stats.normal_flux(2), stats.normal_flux(3));

fprintf('Export:\n');
fprintf('  xyzw file:        %s\n', xyzwFile);
fprintf('  package folder:   %s\n', packageDir);
fprintf('  go3 file:         %s\n', payload_package.surface_file);
fprintf('  gidmsh file:      %s\n', payload_package.scaffold_file);
fprintf('  edges file:       %s\n', payload_package.edges_file);
fprintf('  exported nodes:   %d\n', payload_xyzw.total_nodes);
fprintf('  go3 nodes:        %d\n', payload_package.surface.total_nodes);
fprintf('  scaffold tris:    %d\n', payload_package.scaffold.num_triangles);
fprintf('  CAD edges:        %d\n', payload_package.edges.num_cad_edges);
fprintf('  edge panels:      %d\n\n', payload_package.edges.num_panels);

%% Plots
figure('Name', 'Gmsh scaffold');
surfsmooth3d.stepmesher.plot_scaffold(mesh);
title('Gmsh scaffold');

figure('Name', 'Projected RV nodes');
surfsmooth3d.stepmesher.plot_projected_nodes(mesh, 'filled');
title('Projected RV nodes');

figure('Name', 'Projection distance');
surfsmooth3d.stepmesher.plot_projection_error(mesh);
title('Projection distance');

figure('Name', 'Scaffold minimum angle');
surfsmooth3d.stepmesher.plot_scaffold_quality(mesh, 'min_angle');

figure('Name', 'Scaffold edge ratio');
surfsmooth3d.stepmesher.plot_scaffold_quality(mesh, 'edge_ratio');

figure('Name', 'surfer surface');
surfsmooth3d.stepmesher.plot_surface(mesh);
title('surfer surface');

%% Local helpers
function name = get_filename(pathname)
[~, name, ext] = fileparts(pathname);
name = [name, ext];
end

function stem = erase_file_extension(filename)
[~, stem] = fileparts(filename);
end

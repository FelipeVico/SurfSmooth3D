%% STEP Mesher MATLAB Convergence Driver
% MATLAB analog of the Python step_mesher.ipynb convergence workflow.
%
% This script runs a small parameter study over orders, mesh fractions, and
% refinement levels using the C++/MEX mesher and SurfSmooth3D MATLAB tools.

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
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_geometry_1.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_geometry_2.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_geometry_manas.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'UV_curved_bezier.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_curved_uv.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'cylinder_bezier_v2.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'intersection_two_balls.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'piecewise_bezier_cap.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'spheres_intersect.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'torus_ellipse.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'Multiscale_step.step');


ORDERS = [4, 6];
%MESH_FRACTIONS = [0.1];
MESH_FRACTIONS = [0.05];

REFINEMENT_LEVELS = [0, 1, 2];


EDGE_GLL_ORDER = max(ORDERS);





% Choose one scaffold/adaptivity launch preset for the whole study.
% To test curvature-based Gmsh sizing, change this to 'curvature_balanced'.
% Other useful choices are 'curvature_coarse',
% 'curvature_balanced_local', 'curvature_balanced_local_cadedge',
% 'curvature_balanced_wide', 'curvature_fine', and 'curvature_guarded'.
% The full menu is printed when the driver starts.
STEP1_MESH_PROFILE_CHOICE = 'rigid';
% STEP1_MESH_PROFILE_CHOICE = 'curvature_balanced';
% STEP1_MESH_PROFILE_CHOICE = 'curvature_balanced_local';
% STEP1_MESH_PROFILE_CHOICE = 'curvature_balanced_local_cadedge';
% STEP1_MESH_PROFILE_CHOICE = 'curvature_balanced_wide';
% STEP1_MESH_PROFILE_CHOICE = 'curvature_fine';
% STEP1_MESH_PROFILE_CHOICE = 'curvature_guarded';

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

CURVATURE_POINTS_OVERRIDE = 35;
% HMIN_MODE_OVERRIDE = 'cad_edge';
HMIN_FRACTION_OVERRIDE = 0.01;
% HMIN_EDGE_SCALE_OVERRIDE = 0.25;
% CAD_EDGE_TOL_FRACTION_OVERRIDE = 1e-6;
HMAX_FRACTION_OVERRIDE = 0.12;
MESH_SIZE_EXTEND_FROM_BOUNDARY_OVERRIDE = false;

RUN_PLOTS = true;
SAVE_FIGURES = true;
SAVE_STUDY_CSV = true;
EXPORT_XYZW = false;
EXPORT_GO3 = true;
EXPORT_GIDMSH = true;
EXPORT_GEOM_PACKAGE = true;
KEEP_GOING_AFTER_FAILURE = true;

%% Output folder
stepName = erase_file_extension(get_filename(stepFile));
studyDir = fullfile(repoRoot, 'outputs', 'step', ['convergence_', stepName]);
if ~exist(studyDir, 'dir')
    mkdir(studyDir);
end
figureDir = fullfile(studyDir, 'figures');
if SAVE_FIGURES && ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

fprintf('STEP file:      %s\n', stepFile);
fprintf('Study folder:   %s\n', studyDir);
fprintf('Orders:         %s\n', mat2str(ORDERS));
fprintf('Mesh fractions: %s\n', mat2str(MESH_FRACTIONS));
fprintf('Refinements:    %s\n\n', mat2str(REFINEMENT_LEVELS));

fprintf('Available scaffold/adaptivity choices:\n');
disp(surfsmooth3d.stepmesher.mesh_profile_choices());
[~, studyMeshProfileChoice] = surfsmooth3d.stepmesher.apply_mesh_profile_choice( ...
    struct(), STEP1_MESH_PROFILE_CHOICE);
fprintf('Selected study choice: %s\n', studyMeshProfileChoice.name);
fprintf('  description: %s\n', studyMeshProfileChoice.description);
fprintf('  C++ profile: %s\n\n', studyMeshProfileChoice.mesh_profile);

%% Run convergence study
records = struct([]);
failures = struct([]);
irec = 0;
ifail = 0;

for imf = 1:numel(MESH_FRACTIONS)
    meshFraction = MESH_FRACTIONS(imf);

    fprintf('\n%s\n', repmat('=', 1, 72));
    fprintf('MESH_FRACTION = %.8g\n', meshFraction);
    fprintf('%s\n', repmat('=', 1, 72));

    for iord = 1:numel(ORDERS)
        order = ORDERS(iord);

        for irl = 1:numel(REFINEMENT_LEVELS)
            refinementLevel = REFINEMENT_LEVELS(irl);
            casePrefix = make_case_prefix(order, meshFraction, refinementLevel);

            fprintf('\n--- Case %s ---\n', casePrefix);

            opts = struct();
            opts.order = order;
            opts.mesh_fraction = meshFraction;
            opts.refinement_level = refinementLevel;
            opts.edge_gll_order = EDGE_GLL_ORDER;
            opts = surfsmooth3d.stepmesher.apply_mesh_profile_choice( ...
                opts, STEP1_MESH_PROFILE_CHOICE);
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

            try
                mesh = surfsmooth3d.stepmesher.mesh_step(stepFile, opts);
                stats = surfsmooth3d.stepmesher.validate_area_flux(mesh);

                if EXPORT_XYZW
                    xyzwFile = fullfile(studyDir, [casePrefix, '_xyzw.txt']);
                    payload = surfsmooth3d.stepmesher.export_xyzw(mesh, xyzwFile);
                    exportedNodes = payload.total_nodes;
                else
                    xyzwFile = "";
                    exportedNodes = size(mesh.xyz, 2);
                end

                if EXPORT_GEOM_PACKAGE
                    packageDir = fullfile(studyDir, casePrefix);
                    packagePayload = surfsmooth3d.stepmesher.export_geom_package(mesh, ...
                        packageDir, casePrefix);
                    go3File = string(packagePayload.surface_file);
                    gidmshFile = string(packagePayload.scaffold_file);
                    edgesFile = string(packagePayload.edges_file);
                    exportedNodes = packagePayload.surface.total_nodes;
                elseif EXPORT_GO3 || EXPORT_GIDMSH
                    packageDir = "";
                    go3File = "";
                    gidmshFile = "";
                    edgesFile = "";
                    if EXPORT_GO3
                        go3File = fullfile(studyDir, [casePrefix, '.go3']);
                        go3Payload = surfsmooth3d.stepmesher.export_go3(mesh, go3File);
                        exportedNodes = go3Payload.total_nodes;
                    end
                    if EXPORT_GIDMSH
                        gidmshFile = fullfile(studyDir, ...
                            [casePrefix, '.gidmsh']);
                        surfsmooth3d.stepmesher.export_gidmsh(mesh, gidmshFile);
                    end
                else
                    packageDir = "";
                    gidmshFile = "";
                    go3File = "";
                    edgesFile = "";
                end

                irec = irec + 1;
                records(irec).prefix = string(casePrefix); %#ok<SAGROW>
                records(irec).step_file = string(stepFile);
                records(irec).order = order;
                records(irec).mesh_fraction = meshFraction;
                records(irec).refinement_level = refinementLevel;
                records(irec).edge_gll_order = EDGE_GLL_ORDER;
                records(irec).mesh_profile_choice = ...
                    string(STEP1_MESH_PROFILE_CHOICE);
                records(irec).mesh_profile = string(opts.mesh_profile);
                records(irec).curvature_points = opts.curvature_points;
                records(irec).hmin_mode = string(opts.hmin_mode);
                records(irec).hmin_fraction = opts.hmin_fraction;
                records(irec).hmin_edge_scale = opts.hmin_edge_scale;
                records(irec).cad_edge_tol_fraction = ...
                    opts.cad_edge_tol_fraction;
                records(irec).hmax_fraction = opts.hmax_fraction;
                records(irec).hmin_effective = mesh.meta.hmin_effective;
                records(irec).hmax_effective = mesh.meta.hmax_effective;
                records(irec).min_raw_cad_edge_length = ...
                    mesh.meta.min_raw_cad_edge_length;
                records(irec).min_meaningful_cad_edge_length = ...
                    mesh.meta.min_meaningful_cad_edge_length;
                records(irec).cad_edge_count = mesh.meta.cad_edge_count;
                records(irec).meaningful_cad_edge_count = ...
                    mesh.meta.meaningful_cad_edge_count;
                records(irec).mesh_size_extend_from_boundary = ...
                    opts.mesh_size_extend_from_boundary;
                records(irec).gmsh_algorithm = string(opts.gmsh_algorithm);
                records(irec).optimize = opts.optimize;
                records(irec).min_angle_deg = opts.min_angle_deg;
                records(irec).max_edge_ratio = opts.max_edge_ratio;
                records(irec).enforce_quality = opts.enforce_quality;
                records(irec).min_scaffold_angle_deg = ...
                    mesh.meta.min_scaffold_angle_deg;
                records(irec).max_scaffold_edge_ratio = ...
                    mesh.meta.max_scaffold_edge_ratio;
                records(irec).bad_min_angle_count = ...
                    mesh.meta.bad_min_angle_count;
                records(irec).bad_edge_ratio_count = ...
                    mesh.meta.bad_edge_ratio_count;
                records(irec).n_tri = mesh.npatches;
                records(irec).nodes_per_tri = mesh.nodes_per_patch;
                records(irec).total_nodes = size(mesh.xyz, 2);
                records(irec).bbox_diag = mesh.meta.bbox_diag;
                records(irec).mesh_size = mesh.meta.mesh_size;
                records(irec).h_est = mesh.meta.mesh_size / (2^refinementLevel);
                records(irec).cad_area = stats.cad_area;
                records(irec).numerical_area = stats.area;
                records(irec).abs_area_error = abs(stats.area - stats.cad_area);
                records(irec).rel_area_error = stats.relative_area_error;
                records(irec).normal_flux_x = stats.normal_flux(1);
                records(irec).normal_flux_y = stats.normal_flux(2);
                records(irec).normal_flux_z = stats.normal_flux(3);
                records(irec).normal_flux_norm = norm(stats.normal_flux);
                records(irec).proj_max = mesh.meta.max_projection_distance;
                records(irec).proj_mean = mesh.meta.mean_projection_distance;
                records(irec).projection_failures = mesh.meta.projection_failures;
                records(irec).xyzw_file = string(xyzwFile);
                records(irec).package_dir = string(packageDir);
                records(irec).go3_file = string(go3File);
                records(irec).gidmsh_file = string(gidmshFile);
                records(irec).edges_file = string(edgesFile);
                records(irec).per_face_tags = {stats.per_face_tags(:).'};
                records(irec).per_face_areas = {stats.per_face_areas(:).'};
                records(irec).per_coarse_triangle_ids = {stats.per_coarse_triangle_ids(:).'};
                records(irec).per_coarse_triangle_areas = {stats.per_coarse_triangle_areas(:).'};
                records(irec).scaffold_nodes = {mesh.scaffold_nodes};
                records(irec).scaffold_triangles = {mesh.scaffold_triangles};

                fprintf(['n_tri=%d  nodes=%d  rel_area_err=%.6e  ', ...
                    'flux_norm=%.6e  proj_max=%.6e  min_ang=%.3g  ', ...
                    'max_ratio=%.3g\n'], ...
                    records(irec).n_tri, records(irec).total_nodes, ...
                    records(irec).rel_area_error, records(irec).normal_flux_norm, ...
                    records(irec).proj_max, records(irec).min_scaffold_angle_deg, ...
                    records(irec).max_scaffold_edge_ratio);
                fprintf('exported/available nodes: %d\n', exportedNodes);

            catch ME
                ifail = ifail + 1;
                failures(ifail).prefix = string(casePrefix); %#ok<SAGROW>
                failures(ifail).step_file = string(stepFile);
                failures(ifail).order = order;
                failures(ifail).mesh_fraction = meshFraction;
                failures(ifail).refinement_level = refinementLevel;
                failures(ifail).message = string(ME.message);

                fprintf(2, 'FAILED: %s\n', ME.message);
                if ~KEEP_GOING_AFTER_FAILURE
                    rethrow(ME);
                end
            end
        end
    end
end

%% Tables and exports
if isempty(records)
    error('stepmesher:convergence', 'No successful convergence cases.');
end

dfStudy = struct2table(records);
dfStudy = add_global_self_convergence(dfStudy);
dfPerFace = make_per_face_table(dfStudy);
dfPerFace = add_per_face_self_convergence(dfPerFace);
dfFailures = struct2table_or_empty(failures);

disp(dfStudy);
if ~isempty(dfFailures)
    fprintf('\nFailures:\n');
    disp(dfFailures);
end

if SAVE_STUDY_CSV
    studyCsv = fullfile(studyDir, [stepName, '_convergence_study.csv']);
    dfStudyCsv = remove_cell_columns(dfStudy);
    writetable(dfStudyCsv, studyCsv);
    fprintf('Saved study CSV: %s\n', studyCsv);

    if ~isempty(dfPerFace)
        perFaceCsv = fullfile(studyDir, [stepName, '_per_face_convergence.csv']);
        writetable(dfPerFace, perFaceCsv);
        fprintf('Saved per-face CSV: %s\n', perFaceCsv);
    end

    if ~isempty(dfFailures)
        failureCsv = fullfile(studyDir, [stepName, '_convergence_failures.csv']);
        writetable(dfFailures, failureCsv);
        fprintf('Saved failure CSV: %s\n', failureCsv);
    end
end

%% Plots
if RUN_PLOTS
    setup_plot_style();

    plot_convergence_vs_x(dfStudy, 'rel_area_error', 'h_est', ...
        'estimated h', 'Relative area error vs CAD/OCC', ...
        'Relative area error vs CAD/OCC', 'area_rel_occ_vs_h', ...
        figureDir, SAVE_FIGURES);
    plot_convergence_vs_x(dfStudy, 'rel_area_error', 'n_tri', ...
        'Number of sampled triangles', 'Relative area error vs CAD/OCC', ...
        'Relative area error vs CAD/OCC', 'area_rel_occ_vs_ntri', ...
        figureDir, SAVE_FIGURES);

    plot_convergence_vs_x(dfStudy, 'rel_area_error_vs_finest_global', 'h_est', ...
        'estimated h', '|A_h - A_{ref}| / |A_{ref}|', ...
        'Area self-convergence vs finest run', ...
        'area_self_consistency_global_vs_h', figureDir, SAVE_FIGURES);

    plot_convergence_vs_x(dfStudy, 'normal_flux_norm', 'h_est', ...
        'estimated h', 'Normal flux norm', ...
        'Normal flux norm', 'normal_flux_norm_vs_h', ...
        figureDir, SAVE_FIGURES);
    plot_convergence_vs_x(dfStudy, 'proj_max', 'h_est', ...
        'estimated h', 'Maximum projection distance', ...
        'Maximum projection distance', 'proj_max_vs_h', ...
        figureDir, SAVE_FIGURES);

    plot_error_vs_order(dfStudy, 'rel_area_error', ...
        'Relative area error vs CAD/OCC', ...
        'Area error vs polynomial order', 'area_rel_occ_vs_order', ...
        figureDir, SAVE_FIGURES);
    plot_error_vs_order(dfStudy, 'rel_area_error_vs_finest_global', ...
        '|A_h - A_{ref}| / |A_{ref}|', ...
        'Area self-convergence vs polynomial order', ...
        'area_self_consistency_global_vs_order', figureDir, SAVE_FIGURES);

    if ~isempty(dfPerFace)
        plot_per_face_self_convergence(dfPerFace, figureDir, SAVE_FIGURES);
        plot_per_face_error_vs_order(dfPerFace, figureDir, SAVE_FIGURES);
    end
    plot_coarse_triangle_self_convergence(dfStudy, figureDir, SAVE_FIGURES);
end

%% Local helpers
function prefix = make_case_prefix(order, meshFraction, refinementLevel)
mfToken = strrep(sprintf('%.5f', meshFraction), '.', 'p');
prefix = sprintf('p%d_mf%s_rl%d', order, mfToken, refinementLevel);
end

function name = get_filename(pathname)
[~, name, ext] = fileparts(pathname);
name = [name, ext];
end

function stem = erase_file_extension(filename)
[~, stem] = fileparts(filename);
end

function tbl = struct2table_or_empty(s)
if isempty(s)
    tbl = table();
else
    tbl = struct2table(s);
end
end

function dfStudy = add_global_self_convergence(dfStudy)
refOrder = max(dfStudy.order);
refCandidates = dfStudy(dfStudy.order == refOrder, :);
[~, localIdx] = min(refCandidates.h_est);
refArea = refCandidates.numerical_area(localIdx);

dfStudy.area_ref_global = repmat(refArea, height(dfStudy), 1);
dfStudy.abs_area_error_vs_finest_global = abs(dfStudy.numerical_area - refArea);
if abs(refArea) > 0
    dfStudy.rel_area_error_vs_finest_global = ...
        dfStudy.abs_area_error_vs_finest_global ./ abs(refArea);
else
    dfStudy.rel_area_error_vs_finest_global = NaN(height(dfStudy), 1);
end
end

function dfPerFace = make_per_face_table(dfStudy)
records = struct([]);
irec = 0;
for i = 1:height(dfStudy)
    tags = get_cell_or_value(dfStudy.per_face_tags, i);
    areas = get_cell_or_value(dfStudy.per_face_areas, i);
    tags = double(tags(:));
    areas = double(areas(:));
    for j = 1:numel(tags)
        irec = irec + 1;
        records(irec).prefix = dfStudy.prefix(i); %#ok<AGROW>
        records(irec).order = dfStudy.order(i);
        records(irec).mesh_fraction = dfStudy.mesh_fraction(i);
        records(irec).refinement_level = dfStudy.refinement_level(i);
        records(irec).h_est = dfStudy.h_est(i);
        records(irec).n_tri = dfStudy.n_tri(i);
        records(irec).face_tag = tags(j);
        records(irec).spectral_area = areas(j);
    end
end

if isempty(records)
    dfPerFace = table();
else
    dfPerFace = struct2table(records);
end
end

function value = get_cell_or_value(column, idx)
if iscell(column)
    value = column{idx};
else
    value = column(idx, :);
end
end

function dfPerFace = add_per_face_self_convergence(dfPerFace)
if isempty(dfPerFace)
    return
end

dfPerFace.spectral_area_ref = NaN(height(dfPerFace), 1);
dfPerFace.self_consistency_abs_error = NaN(height(dfPerFace), 1);
dfPerFace.self_consistency_rel_error = NaN(height(dfPerFace), 1);

faceTags = unique(dfPerFace.face_tag).';
for faceTag = faceTags
    faceMask = dfPerFace.face_tag == faceTag;
    sub = dfPerFace(faceMask, :);
    refOrder = max(sub.order);
    refCandidates = sub(sub.order == refOrder, :);
    [~, localIdx] = min(refCandidates.h_est);
    refArea = refCandidates.spectral_area(localIdx);

    idx = find(faceMask);
    dfPerFace.spectral_area_ref(idx) = refArea;
    dfPerFace.self_consistency_abs_error(idx) = ...
        abs(dfPerFace.spectral_area(idx) - refArea);
    if abs(refArea) > 0
        dfPerFace.self_consistency_rel_error(idx) = ...
            dfPerFace.self_consistency_abs_error(idx) ./ abs(refArea);
    end
end
end

function tbl = remove_cell_columns(tbl)
names = tbl.Properties.VariableNames;
drop = false(size(names));
for k = 1:numel(names)
    drop(k) = iscell(tbl.(names{k}));
end
tbl = removevars(tbl, names(drop));
end

function setup_plot_style()
set(groot, 'defaultAxesFontName', 'Times');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultLineLineWidth', 1.6);
set(groot, 'defaultLineMarkerSize', 6);
end

function plot_convergence_vs_x(dfStudy, yField, xField, xLabelText, yLabelText, ...
    plotTitle, fileSuffix, figureDir, saveFigures)
fig = figure('Name', plotTitle);
hold on
orders = unique(dfStudy.order).';
for order = orders
    mask = dfStudy.order == order;
    sub = sortrows(dfStudy(mask, :), xField);
    x = sub.(xField);
    y = sub.(yField);
    valid = x > 0 & y > 0 & isfinite(x) & isfinite(y);
    if any(valid)
        loglog(x(valid), y(valid), '-o', ...
        'DisplayName', sprintf('p=%d', order), ...
        'LineWidth', 1.5, 'MarkerSize', 6);
    end
end
if strcmp(xField, 'h_est')
    set(gca, 'XDir', 'reverse');
end
grid on
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel(xLabelText);
ylabel(yLabelText);
title(plotTitle);
legend('Location', 'best');
hold off
save_figure(fig, figureDir, fileSuffix, saveFigures);
end

function plot_error_vs_order(dfStudy, yField, yLabelText, plotTitle, fileSuffix, ...
    figureDir, saveFigures)
fig = figure('Name', plotTitle);
hold on
refinementLevels = unique(dfStudy.refinement_level).';
for rl = refinementLevels
    mask = dfStudy.refinement_level == rl;
    sub = sortrows(dfStudy(mask, :), {'mesh_fraction', 'order'});
    x = sub.order;
    y = sub.(yField);
    valid = y > 0 & isfinite(y);
    if any(valid)
        semilogy(x(valid), y(valid), '-o', ...
            'DisplayName', sprintf('rl=%d', rl), ...
            'LineWidth', 1.5, 'MarkerSize', 6);
    end
end
grid on
set(gca, 'YScale', 'log');
xticks(unique(dfStudy.order));
xlabel('Polynomial order p');
ylabel(yLabelText);
title(plotTitle);
legend('Location', 'best');
hold off
save_figure(fig, figureDir, fileSuffix, saveFigures);
end

function plot_per_face_self_convergence(dfPerFace, figureDir, saveFigures)
faceTags = unique(dfPerFace.face_tag).';
nfaces = numel(faceTags);
ncols = min(3, nfaces);
nrows = ceil(nfaces / ncols);

fig = figure('Name', 'Per-face area self-convergence vs h');
tiledlayout(nrows, ncols, 'TileSpacing', 'compact');
for faceTag = faceTags
    nexttile
    hold on
    subFace = dfPerFace(dfPerFace.face_tag == faceTag, :);
    orders = unique(subFace.order).';
    for order = orders
        sub = sortrows(subFace(subFace.order == order, :), ...
            {'mesh_fraction', 'refinement_level'});
        err = sub.self_consistency_abs_error;
        valid = err > 0 & isfinite(err);
        if any(valid)
            loglog(sub.h_est(valid), err(valid), '-o', ...
                'DisplayName', sprintf('p=%d', order), ...
                'LineWidth', 1.2, 'MarkerSize', 5);
        end
    end
    set(gca, 'XDir', 'reverse');
    grid on
    title(sprintf('face %d', faceTag));
    xlabel('estimated h');
    ylabel('|A_f - A_{f,ref}|');
    legend('Location', 'best');
    hold off
end
save_figure(fig, figureDir, 'perface_self_consistency_vs_h', saveFigures);
end

function plot_per_face_error_vs_order(dfPerFace, figureDir, saveFigures)
faceTags = unique(dfPerFace.face_tag).';
nfaces = numel(faceTags);
ncols = min(3, nfaces);
nrows = ceil(nfaces / ncols);

fig = figure('Name', 'Per-face area self-convergence vs p');
tiledlayout(nrows, ncols, 'TileSpacing', 'compact');
for faceTag = faceTags
    nexttile
    hold on
    subFace = dfPerFace(dfPerFace.face_tag == faceTag, :);
    refinementLevels = unique(subFace.refinement_level).';
    for rl = refinementLevels
        sub = sortrows(subFace(subFace.refinement_level == rl, :), ...
            {'mesh_fraction', 'order'});
        err = sub.self_consistency_abs_error;
        valid = err > 0 & isfinite(err);
        if any(valid)
            semilogy(sub.order(valid), err(valid), '-o', ...
                'DisplayName', sprintf('rl=%d', rl), ...
                'LineWidth', 1.2, 'MarkerSize', 5);
        end
    end
    grid on
    set(gca, 'YScale', 'log');
    xticks(unique(dfPerFace.order));
    title(sprintf('face %d', faceTag));
    xlabel('Polynomial order p');
    ylabel('|A_f - A_{f,ref}|');
    legend('Location', 'best');
    hold off
end
save_figure(fig, figureDir, 'perface_self_consistency_vs_order', saveFigures);
end

function plot_coarse_triangle_self_convergence(dfStudy, figureDir, saveFigures)
needed = {'per_coarse_triangle_ids', 'per_coarse_triangle_areas', ...
    'scaffold_nodes', 'scaffold_triangles'};
if ~all(ismember(needed, dfStudy.Properties.VariableNames))
    return
end

dfCoarse = make_coarse_triangle_table(dfStudy);
if isempty(dfCoarse)
    return
end

highestOrder = max(dfCoarse.order);
dfCoarse = dfCoarse(dfCoarse.order == highestOrder, :);
meshFractions = unique(dfCoarse.mesh_fraction).';

for meshFraction = meshFractions
    baseMask = abs(dfCoarse.mesh_fraction - meshFraction) <= ...
        1e-12 * max(1, abs(meshFraction));
    subBase = dfCoarse(baseMask, :);
    if isempty(subBase)
        continue
    end

    subBase = add_coarse_triangle_reference(subBase);
    nonRef = subBase(subBase.h_est > subBase.h_ref .* (1.0 + 1e-14), :);
    if isempty(nonRef)
        continue
    end

    triIds = unique(nonRef.coarse_tri_id).';
    bestErr = NaN(numel(triIds), 1);
    for k = 1:numel(triIds)
        err = nonRef.self_consistency_abs_error(nonRef.coarse_tri_id == triIds(k));
        bestErr(k) = min(err, [], 'omitnan');
    end

    plotRows = dfStudy(abs(dfStudy.mesh_fraction - meshFraction) <= ...
        1e-12 * max(1, abs(meshFraction)) & dfStudy.order == highestOrder, :);
    if isempty(plotRows)
        continue
    end
    [~, refRowIdx] = min(plotRows.h_est);
    plotRow = plotRows(refRowIdx, :);

    nodes = get_cell_or_value(plotRow.scaffold_nodes, 1);
    tris = get_cell_or_value(plotRow.scaffold_triangles, 1);
    values = NaN(size(tris, 2), 1);
    for k = 1:numel(triIds)
        triId = triIds(k);
        if triId >= 1 && triId <= numel(values)
            values(triId) = bestErr(k);
        end
    end

    meshForPlot = struct();
    meshForPlot.scaffold_nodes = nodes;
    meshForPlot.scaffold_triangles = tris;

    mfToken = strrep(sprintf('%.5f', meshFraction), '.', 'p');
    figTitle = sprintf(['Coarse-triangle best self-consistency ', ...
        'error, p=%d, mesh fraction %.5g'], highestOrder, meshFraction);
    fig = figure('Name', figTitle);
    ax = axes(fig);
    surfsmooth3d.stepmesher.plot_coarse_triangle_error(meshForPlot, values, ...
        'Parent', ax, ...
        'Title', figTitle, ...
        'ColorbarLabel', 'min |A_T - A_{T,ref}|');

    save_figure(fig, figureDir, ...
        ['coarse_triangle_self_consistency_mf', mfToken], saveFigures);
end
end

function dfCoarse = make_coarse_triangle_table(dfStudy)
records = struct([]);
irec = 0;
for i = 1:height(dfStudy)
    triIds = get_cell_or_value(dfStudy.per_coarse_triangle_ids, i);
    areas = get_cell_or_value(dfStudy.per_coarse_triangle_areas, i);
    triIds = double(triIds(:));
    areas = double(areas(:));
    if numel(triIds) ~= numel(areas)
        continue
    end

    for j = 1:numel(triIds)
        irec = irec + 1;
        records(irec).prefix = dfStudy.prefix(i); %#ok<AGROW>
        records(irec).order = dfStudy.order(i);
        records(irec).mesh_fraction = dfStudy.mesh_fraction(i);
        records(irec).refinement_level = dfStudy.refinement_level(i);
        records(irec).h_est = dfStudy.h_est(i);
        records(irec).n_tri = dfStudy.n_tri(i);
        records(irec).coarse_tri_id = triIds(j);
        records(irec).spectral_area = areas(j);
    end
end

if isempty(records)
    dfCoarse = table();
else
    dfCoarse = struct2table(records);
end
end

function tbl = add_coarse_triangle_reference(tbl)
tbl.spectral_area_ref = NaN(height(tbl), 1);
tbl.h_ref = NaN(height(tbl), 1);
tbl.self_consistency_abs_error = NaN(height(tbl), 1);
tbl.self_consistency_rel_error = NaN(height(tbl), 1);

triIds = unique(tbl.coarse_tri_id).';
for triId = triIds
    mask = tbl.coarse_tri_id == triId;
    sub = tbl(mask, :);
    [~, localIdx] = min(sub.h_est);
    refArea = sub.spectral_area(localIdx);
    hRef = sub.h_est(localIdx);

    idx = find(mask);
    tbl.spectral_area_ref(idx) = refArea;
    tbl.h_ref(idx) = hRef;
    tbl.self_consistency_abs_error(idx) = abs(tbl.spectral_area(idx) - refArea);
    if abs(refArea) > 0
        tbl.self_consistency_rel_error(idx) = ...
            tbl.self_consistency_abs_error(idx) ./ abs(refArea);
    end
end
end

function save_figure(fig, figureDir, fileSuffix, saveFigures)
if ~saveFigures
    return
end
pngFile = fullfile(figureDir, [fileSuffix, '.png']);
pdfFile = fullfile(figureDir, [fileSuffix, '.pdf']);
exportgraphics(fig, pngFile, 'Resolution', 200);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
fprintf('Saved figure: %s\n', pngFile);
fprintf('Saved figure: %s\n', pdfFile);
end

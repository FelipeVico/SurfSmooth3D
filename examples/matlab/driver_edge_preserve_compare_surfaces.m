%DRIVER_EDGE_PRESERVE_COMPARE_SURFACES Compare CAD skeleton and smoothed CAD surface.

driverFile = mfilename('fullpath');
examplesDir = fileparts(driverFile);
repoRoot = fileparts(fileparts(examplesDir));
matlabDir = fullfile(repoRoot, 'matlab');
addpath(matlabDir);
setup_surfsmooth3d;
packageRoot = surfsmooth3d.edgepreserve.default_input_root(repoRoot);

%% User parameters
% Choose one package folder by uncommenting the desired line.
% packageDir = fullfile(packageRoot, ...
%    'convergence_test_geometry_manas', 'p8_mf0p10000_rl0');
% packageDir = fullfile(packageRoot, ...
%     'convergence_test_geometry_manas', 'p8_mf0p10000_rl1');
% packageDir = fullfile(packageRoot, ...
%     'convergence_test_geometry_manas', 'p8_mf0p10000_rl2');
% packageDir = fullfile(packageRoot, ...
%     'convergence_test_geometry_manas', 'p8_mf0p10000_rl3');
% packageDir = fullfile(packageRoot, ...
%     'convergence_test_geometry_2', 'p8_mf0p10000_rl0');

% packageDir = fullfile(packageRoot, ...
%    'convergence_test_geometry_1', 'p8_mf0p10000_rl0');

% packageDir = fullfile(packageRoot, ...
%     'convergence_UV_curved_bezier', 'p8_mf0p10000_rl0');
packageDir = fullfile(repoRoot, 'tests', 'fixtures', 'cad', ...
    'convergence_intersection_two_balls', 'p4_mf0p10000_rl0');

% packageDir = fullfile(packageRoot, ...
%     'convergence_test_curved_uv', 'p8_mf0p05000_rl0');


rlam = 2;
useTwoStageSmoother = true;
makePlots = usejava('desktop');
allowInvalidGo3Save = false;
hideSurferPatchEdges = true;
surfacePlotRefine = 20;
surfacePlotArgs = {'FaceLighting', 'gouraud', ...
    'SpecularStrength', 0.05, ...
    'DiffuseStrength', 0.65, ...
    'AmbientStrength', 0.45};
if hideSurferPatchEdges
    surfacePlotArgs = [{'EdgeColor', 'none', 'LineStyle', 'none'}, ...
        surfacePlotArgs];
end

% Edge-preserving blend controls. In the blend convention below, beta is
% the smooth-surface weight: beta = 0 means CAD, beta = 1 means smooth.
selectedEdgeIds = 2;
c_eps = 0.5;
c_rad = 2.7;
c_tau = 1.05;
sigmaFloorRel = 1e-12;
distanceFloor = 1e-300;

%% Load package inputs
cadFile = fullfile(packageDir, 'surface.go3');
scaffoldFile = fullfile(packageDir, 'scaffold.gidmsh');
edgesFile = fullfile(packageDir, 'edges_gll.txt');

if ~isfile(cadFile)
    error('driver_edge_preserve_compare_surfaces:missingCadFile', ...
        'Could not find CAD skeleton file: %s', cadFile);
end
if ~isfile(scaffoldFile)
    error('driver_edge_preserve_compare_surfaces:missingScaffoldFile', ...
        'Could not find scaffold file: %s', scaffoldFile);
end
if ~isfile(edgesFile)
    error('driver_edge_preserve_compare_surfaces:missingEdgesFile', ...
        'Could not find CAD edge file: %s', edgesFile);
end

S_cad = surfsmooth3d.surfer.load_from_file(cadFile);
edges = load_package_edges(edgesFile);
[scaffoldNodes, scaffoldTris] = read_simple_gidmsh(scaffoldFile);
scaffold = struct();
scaffold.nodes = scaffoldNodes;
scaffold.tris = scaffoldTris;
norder = infer_common_order(S_cad);

tempDir = fullfile(tempdir, 'surfsmooth3d_edge_preserve');
if ~exist(tempDir, 'dir')
    mkdir(tempDir);
end
[~, packageName] = fileparts(packageDir);
tempCadSkeletonTxt = fullfile(tempDir, ...
    sprintf('%s_cad_skeleton_from_go3.txt', packageName));
write_cad_skeleton_txt(S_cad, tempCadSkeletonTxt);

%% Run the CAD-based smoother on the matching scaffold
opts = struct();
opts.filetype = 3;
opts.nrefine = 0;
opts.nquad = norder;
opts.rlam = rlam;
opts.fcad = tempCadSkeletonTxt;
opts.two_stage_smoother = useTwoStageSmoother;

[S_smooth, sigmaCad, gradSigmaCad] = recompute_smooth_state( ...
    scaffoldFile, norder, opts, S_cad, rlam);

%% Report
dispVec = sqrt(sum((S_smooth.r - S_cad.r).^2, 1));
blendParams = struct();
blendParams.selectedEdgeIds = selectedEdgeIds;
blendParams.c_eps = c_eps;
blendParams.c_rad = c_rad;
blendParams.c_tau = c_tau;
blendParams.sigma = sigmaCad;
blendParams.gradSigma = gradSigmaCad;
blendParams.sigmaFloorRel = sigmaFloorRel;
blendParams.geometryDiameter = geometry_diameter(S_cad.r);
blendParams.distanceFloor = distanceFloor;
blendParams.rlam = rlam;
blendParams.useTwoStageSmoother = useTwoStageSmoother;
[S_blend, beta, dsoft] = build_blended_surface_spectral(S_cad, S_smooth, ...
    edges, blendParams);

reportCad = surface_validity_report(S_cad, 'S_cad', []);
reportSmooth = surface_validity_report(S_smooth, 'S_smooth', S_cad);
reportBlend = surface_validity_report(S_blend, 'S_blend', S_cad);
print_surface_validity_report(reportCad);
print_surface_validity_report(reportSmooth);
print_surface_validity_report(reportBlend);
export_bad_surface_state_if_needed(S_cad, S_smooth, S_blend, beta, ...
    dsoft, sigmaCad, blendParams, {reportCad, reportSmooth, reportBlend});

fprintf('driver_edge_preserve_compare_surfaces\n');
fprintf('  packageDir: %s\n', packageDir);
fprintf('  scaffold:   %s\n', scaffoldFile);
fprintf('  CAD go3:    %s\n', cadFile);
fprintf('  edges:      %s\n', edgesFile);
fprintf('  CAD txt:    %s\n', tempCadSkeletonTxt);
fprintf('  norder=%d, npatches=%d, npts=%d, rlam=%g\n', ...
    norder, S_cad.npatches, S_cad.npts, rlam);
fprintf('  two-stage smoother = %d\n', useTwoStageSmoother);
fprintf('  CAD area      = %.16e\n', S_cad.area);
fprintf('  smoothed area = %.16e\n', S_smooth.area);
fprintf('  displacement min/mean/max/std = %.6e %.6e %.6e %.6e\n', ...
    min(dispVec), mean(dispVec), max(dispVec), std(dispVec));
fprintf('  selected edge IDs: %s\n', mat2str(selectedEdgeIds));
fprintf('  c_eps=%g, c_rad=%g, c_tau=%g\n', c_eps, c_rad, c_tau);
fprintf('  sigma min/mean/max = %.6e %.6e %.6e\n', ...
    min(sigmaCad), mean(sigmaCad), max(sigmaCad));
fprintf('  dsoft min/mean/max = %.6e %.6e %.6e\n', ...
    min(dsoft), mean(dsoft), max(dsoft));
fprintf('  beta min/mean/max (0 CAD, 1 smooth) = %.6e %.6e %.6e\n', ...
    min(beta), mean(beta), max(beta));

%% Plots
if makePlots
    figure('Name', 'edge preserve: CAD skeleton');
    plot_surfer_smooth_lit(S_cad, [], surfacePlotArgs, ...
        surfacePlotRefine);
    axis equal;
    view(3);
    camlight headlight;
    colorbar;
    title([packageName ': CAD skeleton'], 'Interpreter', 'none');
    xlabel('x');
    ylabel('y');
    zlabel('z');

    figure('Name', 'edge preserve: smoothed CAD surface');
    plot_surfer_smooth_lit(S_smooth, [], surfacePlotArgs, ...
        surfacePlotRefine);
    axis equal;
    view(3);
    camlight headlight;
    colorbar;
    title([packageName ': smoothed surface from CAD skeleton'], ...
        'Interpreter', 'none');
    xlabel('x');
    ylabel('y');
    zlabel('z');

    figure('Name', 'edge preserve: CAD vs smoothed overlay');
    h1 = plot_surfer_smooth_lit(S_cad, zeros(S_cad.npts, 1), ...
        [surfacePlotArgs, {'FaceColor', [0.10 0.35 0.95], ...
        'FaceAlpha', 0.30}], surfacePlotRefine);
    hold on;
    h2 = plot_surfer_smooth_lit(S_smooth, ones(S_smooth.npts, 1), ...
        [surfacePlotArgs, {'FaceColor', [0.95 0.25 0.10], ...
        'FaceAlpha', 0.30}], surfacePlotRefine);
    axis equal;
    view(3);
    camlight headlight;
    grid on;
    title([packageName ': CAD skeleton and smoothed surface'], ...
        'Interpreter', 'none');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    legend([h1 h2], {'CAD skeleton', 'smoothed'}, 'Location', 'best');

    figure('Name', 'edge preserve: displacement');
    plot_surfer_smooth_lit(S_cad, dispVec, surfacePlotArgs, ...
        surfacePlotRefine);
    axis equal;
    view(3);
    camlight headlight;
    colorbar;
    title([packageName ': |x_{smooth} - x_{CAD}|'], 'Interpreter', 'tex');
    xlabel('x');
    ylabel('y');
    zlabel('z');

    figure('Name', 'edge preserve: beta');
    plot_surfer_smooth_lit(S_cad, beta, surfacePlotArgs, ...
        surfacePlotRefine);
    axis equal;
    view(3);
    camlight headlight;
    colorbar;
    title([packageName ': smooth weight beta (0 CAD, 1 smooth)'], ...
        'Interpreter', 'none');
    xlabel('x');
    ylabel('y');
    zlabel('z');

    figure('Name', 'edge preserve: blended surface');
    plot_surfer_smooth_lit(S_blend, beta, surfacePlotArgs, ...
        surfacePlotRefine);
    hold on;
    plot_selected_edges(gca, edges, selectedEdgeIds, [0 0 0], 2.5);
    axis equal;
    view(3);
    camlight headlight;
    colorbar;
    title([packageName ': beta*smooth + (1-beta)*CAD'], ...
        'Interpreter', 'none');
    xlabel('x');
    ylabel('y');
    zlabel('z');

    launch_interactive_blend_viewer(S_cad, S_smooth, edges, ...
        blendParams, packageName, surfacePlotArgs, surfacePlotRefine, ...
        scaffoldFile, scaffold, norder, opts, repoRoot, ...
        allowInvalidGo3Save);
else
    fprintf('  plots skipped because MATLAB is running without desktop graphics\n');
end

function norder = infer_common_order(S)
    norders = S.norders(:);
    if isempty(norders)
        error('driver_edge_preserve_compare_surfaces:emptyOrders', ...
            'CAD skeleton has no patch orders.');
    end
    if any(norders ~= norders(1))
        error('driver_edge_preserve_compare_surfaces:mixedOrders', ...
            'This driver expects one common order across all patches.');
    end
    norder = norders(1);
end

function write_cad_skeleton_txt(S, outFile)
    [srcvals, ~, ~, ~, ~, wts] = extract_arrays(S);
    data = [srcvals(1:3, :).', srcvals(10:12, :).', wts(:)];

    fid = fopen(outFile, 'w');
    if fid < 0
        error('driver_edge_preserve_compare_surfaces:openFailed', ...
            'Could not open temporary CAD skeleton file for writing: %s', ...
            outFile);
    end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '# x y z nx ny nz w generated from surface.go3\n');
    fprintf(fid, '%.16e %.16e %.16e %.16e %.16e %.16e %.16e\n', data.');
    clear cleanup;
end

function [nodes, tris] = read_simple_gidmsh(fname)
    lines = string(readlines(fname));
    iCoord = find(contains(lines, "Coordinates"), 1, 'first');
    iEndCoord = find(contains(lines, "End Coordinates"), 1, 'first');
    iElem = find(contains(lines, "Elements"), 1, 'first');
    iEndElem = find(contains(lines, "End Elements"), 1, 'first');
    if isempty(iCoord) || isempty(iEndCoord) || isempty(iElem) || ...
            isempty(iEndElem)
        error('driver_edge_preserve_compare_surfaces:badGidmsh', ...
            'Could not parse GiD mesh sections in %s.', fname);
    end

    nodeLines = lines((iCoord + 1):(iEndCoord - 1));
    elemLines = lines((iElem + 1):(iEndElem - 1));
    nodeData = sscanf(join(nodeLines, newline), '%f');
    elemData = sscanf(join(elemLines, newline), '%f');
    if mod(numel(nodeData), 4) ~= 0 || mod(numel(elemData), 4) ~= 0
        error('driver_edge_preserve_compare_surfaces:badGidmshData', ...
            'Unexpected GiD data length in %s.', fname);
    end

    nodeData = reshape(nodeData, 4, []);
    elemData = reshape(elemData, 4, []);
    ids = nodeData(1, :);
    if any(ids ~= 1:numel(ids))
        error('driver_edge_preserve_compare_surfaces:noncontiguousNodes', ...
            'This viewer expects contiguous 1-based GiD node ids.');
    end
    nodes = nodeData(2:4, :);
    tris = elemData(2:4, :);
end

function assert_compatible_surfaces(S_cad, S_smooth)
    if S_cad.npatches ~= S_smooth.npatches
        error('driver_edge_preserve_compare_surfaces:patchMismatch', ...
            ['CAD and smoothed surfaces have different patch counts: ', ...
            '%d vs %d. Use opts.nrefine = 0 and matching package files.'], ...
            S_cad.npatches, S_smooth.npatches);
    end
    if S_cad.npts ~= S_smooth.npts
        error('driver_edge_preserve_compare_surfaces:nodeMismatch', ...
            'CAD and smoothed surfaces have different node counts: %d vs %d.', ...
            S_cad.npts, S_smooth.npts);
    end
    if ~isequal(S_cad.norders(:), S_smooth.norders(:))
        error('driver_edge_preserve_compare_surfaces:orderMismatch', ...
            'CAD and smoothed surfaces have different patch orders.');
    end
    if ~isequal(S_cad.iptype(:), S_smooth.iptype(:))
        error('driver_edge_preserve_compare_surfaces:iptypeMismatch', ...
            'CAD and smoothed surfaces have different patch types.');
    end
end

function report = surface_validity_report(S, label, S_ref)
    if nargin < 3
        S_ref = [];
    end

    thresholds = struct();
    thresholds.normal_cross_dot_min = 0.99;
    thresholds.jac_ratio_fail_min = 1e-3;
    thresholds.jac_ratio_warn_min = 1e-2;

    [srcvals, ~, ~, ~, ~, wts] = extract_arrays(S);
    cr = cross(srcvals(4:6, :).', srcvals(7:9, :).').';
    jac = sqrt(sum(cr.^2, 1)).';
    crUnit = cr ./ max(jac.', realmin);
    normalCrossDot = sum(srcvals(10:12, :) .* crUnit, 1).';

    hasRef = ~isempty(S_ref) && S_ref.npts == S.npts && ...
        S_ref.npatches == S.npatches;
    if hasRef
        [srcRef, ~, ~, ~, ~, ~] = extract_arrays(S_ref);
        crRef = cross(srcRef(4:6, :).', srcRef(7:9, :).').';
        jacRef = sqrt(sum(crRef.^2, 1)).';
        jacRatioToRef = jac ./ max(jacRef, realmin);
        displacementToRef = sqrt(sum((S.r - S_ref.r).^2, 1)).';
        curvatureThreshold = 100 * max(1, max(abs(S_ref.mean_curv(:))));
    else
        jacRatioToRef = nan(S.npts, 1);
        displacementToRef = nan(S.npts, 1);
        curvatureThreshold = inf;
    end
    thresholds.curvature_abs_max = curvatureThreshold;

    finiteNode = all(isfinite(srcvals), 1).' & isfinite(wts(:)) & ...
        isfinite(S.mean_curv(:)) & isfinite(jac) & ...
        isfinite(normalCrossDot);
    if hasRef
        finiteNode = finiteNode & isfinite(jacRatioToRef) & ...
            isfinite(displacementToRef);
    end

    badFinite = ~finiteNode;
    badNormal = normalCrossDot < thresholds.normal_cross_dot_min;
    badJacRatio = hasRef & ...
        (jacRatioToRef < thresholds.jac_ratio_fail_min);
    warnJacRatio = hasRef & ...
        (jacRatioToRef < thresholds.jac_ratio_warn_min) & ...
        ~badJacRatio;
    badCurvature = abs(S.mean_curv(:)) > thresholds.curvature_abs_max;
    badMask = badFinite | badNormal | badJacRatio | badCurvature;

    failureCode = zeros(S.npts, 1);
    failureCode(warnJacRatio) = 1;
    failureCode(badCurvature) = 2;
    failureCode(badJacRatio) = 3;
    failureCode(badNormal) = 4;
    failureCode(badFinite) = 5;

    report = struct();
    report.label = char(label);
    report.thresholds = thresholds;
    report.has_reference = hasRef;
    report.normal_cross_dot = normalCrossDot;
    report.jac = jac;
    report.jac_ratio_to_ref = jacRatioToRef;
    report.displacement_to_ref = displacementToRef;
    report.r = S.r;
    report.patch_id = S.patch_id(:);
    report.wts = wts(:);
    report.mean_curv = S.mean_curv(:);
    report.bad_mask = badMask;
    report.warning_mask = warnJacRatio;
    report.failure_code = failureCode;
    report.is_valid = ~any(badMask);
    report.has_warnings = any(warnJacRatio);
    report.csv_file = fullfile(tempdir,'surfsmooth3d_surface_validity_bad_nodes.csv');

    report.stats = struct();
    report.stats.npts = S.npts;
    report.stats.npatches = S.npatches;
    report.stats.normal_cross_dot_min = min(normalCrossDot);
    report.stats.normal_cross_dot_mean = mean(normalCrossDot);
    report.stats.normal_cross_dot_max = max(normalCrossDot);
    report.stats.jac_min = min(jac);
    report.stats.jac_mean = mean(jac);
    report.stats.jac_max = max(jac);
    report.stats.wts_min = min(wts);
    report.stats.wts_mean = mean(wts);
    report.stats.wts_max = max(wts);
    report.stats.mean_curv_min = min(S.mean_curv);
    report.stats.mean_curv_mean = mean(S.mean_curv);
    report.stats.mean_curv_max = max(S.mean_curv);
    report.stats.mean_curv_std = std(S.mean_curv);
    report.stats.mean_curv_abs_max = max(abs(S.mean_curv));
    report.stats.bad_node_count = nnz(badMask);
    report.stats.warning_node_count = nnz(warnJacRatio);
    report.stats.bad_patch_count = count_patches(S, badMask);
    report.stats.warning_patch_count = count_patches(S, warnJacRatio);
    report.stats.bad_finite_count = nnz(badFinite);
    report.stats.bad_normal_count = nnz(badNormal);
    report.stats.bad_jac_ratio_count = nnz(badJacRatio);
    report.stats.warn_jac_ratio_count = nnz(warnJacRatio);
    report.stats.bad_curvature_count = nnz(badCurvature);
    if hasRef
        report.stats.jac_ratio_min = min(jacRatioToRef);
        report.stats.jac_ratio_mean = mean(jacRatioToRef);
        report.stats.jac_ratio_max = max(jacRatioToRef);
        report.stats.displacement_min = min(displacementToRef);
        report.stats.displacement_mean = mean(displacementToRef);
        report.stats.displacement_max = max(displacementToRef);
        report.stats.displacement_std = std(displacementToRef);
    else
        report.stats.jac_ratio_min = nan;
        report.stats.jac_ratio_mean = nan;
        report.stats.jac_ratio_max = nan;
        report.stats.displacement_min = nan;
        report.stats.displacement_mean = nan;
        report.stats.displacement_max = nan;
        report.stats.displacement_std = nan;
    end
end

function print_surface_validity_report(report)
    status = 'PASS';
    if ~report.is_valid
        status = 'FAIL';
    elseif report.has_warnings
        status = 'WARN';
    end

    fprintf('\nBIE surface validity report: %s [%s]\n', ...
        report.label, status);
    fprintf('  npts=%d, npatches=%d\n', ...
        report.stats.npts, report.stats.npatches);
    fprintf('  n dot normalize(du x dv) min/mean/max = %.6e %.6e %.6e\n', ...
        report.stats.normal_cross_dot_min, ...
        report.stats.normal_cross_dot_mean, ...
        report.stats.normal_cross_dot_max);
    fprintf('  |du x dv| min/mean/max = %.6e %.6e %.6e\n', ...
        report.stats.jac_min, report.stats.jac_mean, ...
        report.stats.jac_max);
    if report.has_reference
        fprintf('  jac/ref min/mean/max = %.6e %.6e %.6e\n', ...
            report.stats.jac_ratio_min, report.stats.jac_ratio_mean, ...
            report.stats.jac_ratio_max);
        fprintf('  displacement to ref min/mean/max/std = %.6e %.6e %.6e %.6e\n', ...
            report.stats.displacement_min, report.stats.displacement_mean, ...
            report.stats.displacement_max, report.stats.displacement_std);
    end
    fprintf('  weights min/mean/max = %.6e %.6e %.6e\n', ...
        report.stats.wts_min, report.stats.wts_mean, ...
        report.stats.wts_max);
    fprintf('  mean curvature min/mean/max/std = %.6e %.6e %.6e %.6e\n', ...
        report.stats.mean_curv_min, report.stats.mean_curv_mean, ...
        report.stats.mean_curv_max, report.stats.mean_curv_std);
    fprintf('  bad nodes/patches = %d/%d, warnings nodes/patches = %d/%d\n', ...
        report.stats.bad_node_count, report.stats.bad_patch_count, ...
        report.stats.warning_node_count, report.stats.warning_patch_count);
    if ~report.is_valid
        fprintf(['  failure counts: nonfinite=%d, normal=%d, ', ...
            'jac_ratio=%d, curvature=%d\n'], ...
            report.stats.bad_finite_count, report.stats.bad_normal_count, ...
            report.stats.bad_jac_ratio_count, ...
            report.stats.bad_curvature_count);
    end
end

function assert_surface_valid_for_bie(report, allowInvalidSave)
    if ~report.is_valid || report.has_warnings
        write_surface_validity_bad_nodes_csv(report);
    end
    if report.is_valid
        return
    end

    msg = sprintf(['Surface "%s" failed BIE .go3 validity checks. ', ...
        'Refusing to save invalid geometry. Bad nodes=%d, bad patches=%d. ', ...
        'Diagnostics CSV: %s'], report.label, ...
        report.stats.bad_node_count, report.stats.bad_patch_count, ...
        report.csv_file);
    if allowInvalidSave
        warning('driver_edge_preserve_compare_surfaces:invalidGo3Save', ...
            '%s\nallowInvalidGo3Save=true, so save will continue.', msg);
    else
        error('driver_edge_preserve_compare_surfaces:invalidGo3Save', ...
            '%s', msg);
    end
end

function write_surface_validity_bad_nodes_csv(report)
    T = surface_validity_bad_nodes_table(report);
    if isempty(T)
        return
    end
    writetable(T, report.csv_file);
    fprintf('  wrote validity diagnostics CSV: %s\n', report.csv_file);
end

function write_surface_validity_reports_csv(reports, csvFile)
    tables = {};
    for i = 1:numel(reports)
        T = surface_validity_bad_nodes_table(reports{i});
        if ~isempty(T)
            tables{end+1} = T; %#ok<AGROW>
        end
    end
    if isempty(tables)
        return
    end
    Tall = vertcat(tables{:});
    writetable(Tall, csvFile);
    fprintf('  wrote combined validity diagnostics CSV: %s\n', csvFile);
end

function T = surface_validity_bad_nodes_table(report)
    mask = report.bad_mask | report.warning_mask;
    if ~any(mask)
        T = table();
        return
    end

    nodeIndex = find(mask);
    nrows = numel(nodeIndex);
    failureLabel = strings(nrows, 1);
    for i = 1:nrows
        failureLabel(i) = failure_code_label(report.failure_code(nodeIndex(i)));
    end

    T = table( ...
        repmat(string(report.label), nrows, 1), ...
        nodeIndex(:), ...
        report.patch_id(nodeIndex), ...
        report.r(1, nodeIndex).', ...
        report.r(2, nodeIndex).', ...
        report.r(3, nodeIndex).', ...
        report.normal_cross_dot(nodeIndex), ...
        report.jac(nodeIndex), ...
        report.jac_ratio_to_ref(nodeIndex), ...
        report.mean_curv(nodeIndex), ...
        report.displacement_to_ref(nodeIndex), ...
        failureLabel, ...
        'VariableNames', {'surface_label', 'node_index', ...
        'patch_index', 'x', 'y', 'z', 'normal_cross_dot', 'jac', ...
        'jac_ratio_to_ref', 'mean_curv', 'displacement_to_ref', ...
        'failure_label'});
end

function txt = failure_code_label(code)
    switch code
        case 5
            txt = "nonfinite";
        case 4
            txt = "normal_cross";
        case 3
            txt = "jac_ratio_fail";
        case 2
            txt = "curvature";
        case 1
            txt = "jac_ratio_warn";
        otherwise
            txt = "ok";
    end
end

function n = count_patches(S, mask)
    if any(mask)
        n = numel(unique(S.patch_id(mask)));
    else
        n = 0;
    end
end

function export_bad_surface_state_if_needed(S_cad, S_smooth, S_blend, ...
        beta, dsoft, sigmaCad, params, reports)
    isBad = false;
    for i = 1:numel(reports)
        isBad = isBad || ~reports{i}.is_valid;
    end
    write_surface_validity_reports_csv(reports, ...
        fullfile(tempdir,'surfsmooth3d_surface_validity_bad_nodes.csv'));
    if ~isBad
        return
    end
    outFile = fullfile(tempdir,'surfsmooth3d_bad_surface_state.mat');
    save(outFile, 'S_cad', 'S_smooth', 'S_blend', 'beta', 'dsoft', ...
        'sigmaCad', 'params', 'reports', '-v7.3');
    fprintf('  wrote bad surface state MAT: %s\n', outFile);
end

function edges = load_package_edges(edgeFile)
    loader = which('surfsmooth3d.stepmesher.load_edges_gll');
    if ~isempty(loader)
        edges = surfsmooth3d.stepmesher.load_edges_gll(edgeFile);
        return
    end

    fid = fopen(edgeFile, 'r');
    if fid < 0
        error('driver_edge_preserve_compare_surfaces:edgesOpenFailed', ...
            'Could not open edge file: %s', edgeFile);
    end
    cleanup = onCleanup(@() fclose(fid));

    marker = next_nonempty_line(fid);
    switch marker
        case 'STEP_MESHER_EDGES_GLL_V2'
            edges = read_edges_v2(fid, edgeFile);
        case 'STEP_MESHER_EDGES_GLL_V1'
            edges = read_edges_v1_as_v2(fid, edgeFile);
        otherwise
            error('driver_edge_preserve_compare_surfaces:badEdgesFile', ...
                'Unsupported edge file marker in %s.', edgeFile);
    end
    clear cleanup;
end

function [S_smooth, sigmaCad, gradSigmaCad] = recompute_smooth_state( ...
        scaffoldFile, norder, opts, S_cad, rlam)
    optsUse = opts;
    optsUse.rlam = rlam;

    Sall = surfsmooth3d.multiscale_mesher(scaffoldFile, norder, optsUse);
    S_smooth = Sall{end};
    assert_compatible_surfaces(S_cad, S_smooth);

    [sigmaCad, gradSigmaCad] = surfsmooth3d.multiscale_mesher_sigma_eval( ...
        scaffoldFile, norder, optsUse, S_cad.r);
    sigmaCad = sigmaCad(:);
    gradSigmaCad = double(gradSigmaCad);
end

function [S_blend, beta, dsoft] = build_blended_surface_spectral( ...
        S_cad, S_smooth, edges, params)
    dsoft = compute_edge_distance(S_cad.r, edges, params);
    beta = beta_from_distance(dsoft, params);

    [srcCad, ~, norders, ~, iptype, ~] = extract_arrays(S_cad);
    [srcSmooth, ~, ~, ~, ~, ~] = extract_arrays(S_smooth);

    betaRow = reshape(beta, 1, []);
    beta3 = repmat(betaRow, 3, 1);
    rBlend = beta3 .* srcSmooth(1:3, :) + ...
        (1 - beta3) .* srcCad(1:3, :);
    srcBlend = srcvals_from_positions(S_cad, rBlend);

    S_blend = surfsmooth3d.surfer(S_cad.npatches, norders, srcBlend, iptype);
end

function srcvals = srcvals_from_positions(S_template, rvals)
    [~, ~, norders, ~, iptype, ~] = extract_arrays(S_template);
    if any(iptype(:) ~= 1)
        error('driver_edge_preserve_compare_surfaces:unsupportedPatchType', ...
            'Spectral blend export currently supports triangular RV patches only.');
    end
    if ~isequal(size(rvals), size(S_template.r))
        error('driver_edge_preserve_compare_surfaces:badBlendPositions', ...
            'Blended position array has the wrong size.');
    end

    srcvals = zeros(12, S_template.npts);
    srcvals(1:3, :) = rvals;

    for ipatch = 1:S_template.npatches
        order = double(norders(ipatch));
        i1 = double(S_template.ixyzs(ipatch));
        i2 = double(S_template.ixyzs(ipatch + 1)) - 1;
        nodeRange = i1:i2;

        rnodes = surfsmooth3d.internal.koorn.rv_nodes(order);
        valsToCoefs = surfsmooth3d.internal.koorn.vals2coefs(order, rnodes);
        [~, dersu, dersv] = surfsmooth3d.internal.koorn.ders(order, rnodes);

        rcoefs = rvals(:, nodeRange) * valsToCoefs.';
        du = rcoefs * dersu;
        dv = rcoefs * dersv;
        n = cross(du, dv, 1);
        nNorm = sqrt(sum(n.^2, 1));
        bad = ~isfinite(nNorm) | nNorm <= eps;
        if any(bad)
            error('driver_edge_preserve_compare_surfaces:degenerateBlendPatch', ...
                'Blended patch %d has degenerate spectral derivatives.', ipatch);
        end
        n = n ./ nNorm;

        srcvals(4:6, nodeRange) = du;
        srcvals(7:9, nodeRange) = dv;
        srcvals(10:12, nodeRange) = n;
    end
end

function beta = beta_from_distance(dsoft, params)
    sigma = local_sigma(params);
    radius = max(params.c_rad * sigma, eps);
    betaTau = max(params.c_tau * sigma, eps);
    beta = 0.5 * erfc((dsoft(:) - radius) ./ (sqrt(2) * betaTau));
    beta = min(max(beta, 0), 1);
end

function dsoft = compute_edge_distance(xyz, edges, params)
    [edgeXyz, edgeW] = selected_edge_quadrature(edges, ...
        params.selectedEdgeIds);
    if isempty(edgeXyz)
        error('driver_edge_preserve_compare_surfaces:noSelectedEdges', ...
            'No valid selected CAD edges were found.');
    end

    sigma = local_sigma(params);
    epsSoftmin = params.c_eps * sigma;
    if any(~isfinite(epsSoftmin)) || any(epsSoftmin <= 0)
        error('driver_edge_preserve_compare_surfaces:badCEps', ...
            'c_eps*sigma must be positive and finite at every node.');
    end

    npts = size(xyz, 2);
    A = zeros(npts, 1);
    chunkSize = 2000;

    for i1 = 1:chunkSize:npts
        i2 = min(npts, i1 + chunkSize - 1);
        xloc = xyz(:, i1:i2).';
        epsLoc = epsSoftmin(i1:i2);
        dx = xloc(:, 1) - edgeXyz(1, :);
        dy = xloc(:, 2) - edgeXyz(2, :);
        dz = xloc(:, 3) - edgeXyz(3, :);
        dist2 = dx.^2 + dy.^2 + dz.^2;
        denom = 2 * epsLoc.^2;
        normalizer = 1 ./ (sqrt(2*pi) * epsLoc);
        A(i1:i2) = normalizer .* ...
            (exp(-dist2 ./ denom) * edgeW(:));
    end

    floorVal = max(params.distanceFloor, realmin);
    dsoft = epsSoftmin .* sqrt(max(0, -2 * log(max(A, floorVal))));
end

function sigma = local_sigma(params)
    sigma = params.sigma(:);
    sigmaFloor = max(params.sigmaFloorRel * params.geometryDiameter, realmin);
    sigma = max(sigma, sigmaFloor);
    if any(~isfinite(sigma)) || any(sigma <= 0)
        error('driver_edge_preserve_compare_surfaces:badSigma', ...
            'The smoother sigma evaluator returned invalid sigma values.');
    end
end

function diam = geometry_diameter(xyz)
    mins = min(xyz, [], 2);
    maxs = max(xyz, [], 2);
    diam = norm(maxs - mins);
    diam = max(diam, 1);
end

function [edgeXyz, edgeW] = selected_edge_quadrature(edges, edgeIds)
    edgeXyz = zeros(3, 0);
    edgeW = zeros(1, 0);
    edgeIds = unique(edgeIds(:).');

    for id = edgeIds
        idx = find([edges.edge_id] == id, 1);
        if isempty(idx)
            continue
        end
        panels = edges(idx).panels;
        for j = 1:numel(panels)
            edgeXyz = [edgeXyz, panels(j).xyz]; %#ok<AGROW>
            edgeW = [edgeW, panels(j).w]; %#ok<AGROW>
        end
    end
end

function h = plot_selected_edges(ax, edges, edgeIds, color, lineWidth)
    h = gobjects(1, 0);
    hold(ax, 'on');
    edgeIds = unique(edgeIds(:).');
    for id = edgeIds
        idx = find([edges.edge_id] == id, 1);
        if isempty(idx)
            continue
        end
        panels = edges(idx).panels;
        for j = 1:numel(panels)
            xyz = panels(j).xyz;
            h(end+1) = line(ax, xyz(1, :), xyz(2, :), xyz(3, :), ...
                'Color', color, 'LineWidth', lineWidth, ...
                'HitTest', 'off'); %#ok<AGROW>
        end
    end
end

function h = plot_flat_scaffold(ax, scaffold, mode)
    hold(ax, 'on');
    faces = double(scaffold.tris.');
    nodes = double(scaffold.nodes);
    switch mode
        case 'solid'
            h = trisurf(faces, nodes(1, :), nodes(2, :), nodes(3, :), ...
                zeros(1, size(nodes, 2)), 'Parent', ax, ...
                'FaceColor', [0.72 0.76 0.80], ...
                'FaceAlpha', 0.78, ...
                'EdgeColor', [0.05 0.05 0.05], ...
                'EdgeAlpha', 0.45, ...
                'LineWidth', 0.6, ...
                'AmbientStrength', 0.55, ...
                'DiffuseStrength', 0.55, ...
                'SpecularStrength', 0.05);
        case 'wireframe'
            h = trisurf(faces, nodes(1, :), nodes(2, :), nodes(3, :), ...
                zeros(1, size(nodes, 2)), 'Parent', ax, ...
                'FaceColor', 'none', ...
                'EdgeColor', [0.02 0.02 0.02], ...
                'EdgeAlpha', 0.28, ...
                'LineWidth', 0.45, ...
                'HitTest', 'off');
        otherwise
            error('driver_edge_preserve_compare_surfaces:badScaffoldMode', ...
                'Unknown scaffold plot mode "%s".', mode);
    end
end

function h = plot_surfer_smooth_lit(S, values, plotArgs, refine)
    if nargin < 2 || isempty(values)
        values = S.mean_curv;
    end
    if nargin < 3 || isempty(plotArgs)
        plotArgs = {};
    end
    if nargin < 4 || isempty(refine)
        refine = 20;
    end

    if any(S.iptype(:) ~= 1)
        h = plot(S, values, plotArgs{:});
        return
    end

    values = expand_surfer_plot_values(S, values);
    fcoefs = S.vals2coefs(values.');
    ncoefs = S.vals2coefs(S.n);
    [uvPlot, triFaces] = local_reference_triangle_mesh(refine);

    npatches = S.npatches;
    nplot = size(uvPlot, 2);
    nfaces = size(triFaces, 1);
    vertices = zeros(npatches*nplot, 3);
    vertexNormals = zeros(npatches*nplot, 3);
    cdata = zeros(npatches*nplot, 1);
    faces = zeros(npatches*nfaces, 3);

    polCache = struct();
    for ipatch = 1:npatches
        order = double(S.norders(ipatch));
        cacheKey = sprintf('o%d', order);
        if ~isfield(polCache, cacheKey)
            polCache.(cacheKey) = surfsmooth3d.internal.koorn.pols(order, uvPlot);
        end
        pols = polCache.(cacheKey);

        i1 = double(S.ixyzs(ipatch));
        i2 = double(S.ixyzs(ipatch + 1)) - 1;
        nodeRange = i1:i2;

        rPlot = S.srccoefs{ipatch}(1:3, :) * pols;
        nPlot = ncoefs(:, nodeRange) * pols;
        nNorm = sqrt(sum(nPlot.^2, 1));
        good = isfinite(nNorm) & nNorm > 0;
        nPlot(:, good) = nPlot(:, good) ./ nNorm(good);
        nPlot(:, ~good) = 0;
        fPlot = fcoefs(nodeRange) * pols;

        vertexRange = (ipatch - 1)*nplot + (1:nplot);
        faceRange = (ipatch - 1)*nfaces + (1:nfaces);
        vertices(vertexRange, :) = rPlot.';
        vertexNormals(vertexRange, :) = nPlot.';
        cdata(vertexRange) = fPlot(:);
        faces(faceRange, :) = triFaces + vertexRange(1) - 1;
    end

    h = patch('Faces', faces, 'Vertices', vertices, ...
        'FaceVertexCData', cdata, ...
        'FaceColor', 'interp', ...
        'EdgeColor', 'none', ...
        'LineStyle', 'none', ...
        'FaceLighting', 'gouraud', ...
        'SpecularStrength', 0.05, ...
        'DiffuseStrength', 0.65, ...
        'AmbientStrength', 0.45, ...
        plotArgs{:});
    set(h, 'VertexNormals', vertexNormals);
    axis equal;
    view(3);
    grid on;
end

function values = expand_surfer_plot_values(S, values)
    values = values(:);
    if numel(values) == S.npts
        return
    end
    if numel(values) == S.npatches
        patchValues = values;
        values = zeros(S.npts, 1);
        for ipatch = 1:S.npatches
            i1 = double(S.ixyzs(ipatch));
            i2 = double(S.ixyzs(ipatch + 1)) - 1;
            values(i1:i2) = patchValues(ipatch);
        end
        return
    end
    error('driver_edge_preserve_compare_surfaces:badPlotValues', ...
        'Plot values must have length S.npts or S.npatches.');
end

function [uv, faces] = local_reference_triangle_mesh(refine)
    refine = max(1, round(refine));
    index = nan(refine + 1, refine + 1);
    uv = zeros(2, (refine + 1)*(refine + 2)/2);
    k = 0;
    for i = 0:refine
        for j = 0:(refine - i)
            k = k + 1;
            index(i + 1, j + 1) = k;
            uv(:, k) = [i; j] / refine;
        end
    end

    faces = zeros(refine*refine, 3);
    nf = 0;
    for i = 0:(refine - 1)
        for j = 0:(refine - 1 - i)
            v00 = index(i + 1, j + 1);
            v10 = index(i + 2, j + 1);
            v01 = index(i + 1, j + 2);
            nf = nf + 1;
            faces(nf, :) = [v00 v10 v01];
            if j <= refine - 2 - i
                v11 = index(i + 2, j + 2);
                nf = nf + 1;
                faces(nf, :) = [v10 v11 v01];
            end
        end
    end
    faces = faces(1:nf, :);
end

function launch_interactive_blend_viewer(S_cad, S_smooth, edges, ...
        params, packageName, surfacePlotArgs, surfacePlotRefine, ...
        scaffoldFile, scaffold, norder, opts, repoRoot, ...
        allowInvalidGo3Save)
    fig = figure('Name', 'edge preserve: interactive blend');
    ax = axes('Parent', fig, 'Units', 'normalized', ...
        'Position', [0.06 0.08 0.66 0.86]);
    panel = uipanel('Parent', fig, 'Title', 'Edge-preserve blend', ...
        'Units', 'normalized', 'Position', [0.75 0.08 0.23 0.86]);

    state = struct();
    state.surfaceHandle = gobjects(0);
    state.edgeHandles = gobjects(0);
    state.scaffoldHandle = gobjects(0);
    state.beta = [];
    state.dsoft = [];
    state.S_smooth = S_smooth;
    state.S_blend = [];
    state.reportSmooth = surface_validity_report(S_smooth, ...
        'S_smooth_live', S_cad);
    state.reportBlend = [];
    state.pendingRlam = params.rlam;
    state.pendingTwoStageSmoother = logical(params.useTwoStageSmoother);
    state.displayMode = 'blend';
    state.axisLimits = geometry_axis_limits(S_cad, state.S_smooth, ...
        scaffold);
    initialParams = params;

    make_controls(panel);
    refresh_plot();

    function make_controls(parent)
        uicontrol(parent, 'Style', 'text', 'String', 'Edge IDs', ...
            'HorizontalAlignment', 'left', 'Units', 'normalized', ...
            'Position', [0.06 0.93 0.88 0.04]);
        state.edgeEdit = uicontrol(parent, 'Style', 'edit', ...
            'String', mat2str(params.selectedEdgeIds), ...
            'HorizontalAlignment', 'left', 'Units', 'normalized', ...
            'Position', [0.06 0.89 0.88 0.045], ...
            'Callback', @(~, ~) apply_from_controls());

        uicontrol(parent, 'Style', 'text', 'String', 'Display mode', ...
            'HorizontalAlignment', 'left', 'Units', 'normalized', ...
            'Position', [0.06 0.835 0.88 0.035]);
        state.displayPopup = uicontrol(parent, 'Style', 'popupmenu', ...
            'String', {'edge-preserve blend', 'fully smooth', ...
            'flat scaffold'}, ...
            'Value', 1, 'Units', 'normalized', ...
            'Position', [0.06 0.805 0.88 0.045], ...
            'Callback', @(~, ~) display_mode_changed());
        state.twoStageCheckbox = uicontrol(parent, 'Style', 'checkbox', ...
            'String', 'two-stage smoother (press Draw/Recompute)', ...
            'Value', double(params.useTwoStageSmoother), ...
            'Units', 'normalized', 'Position', [0.06 0.755 0.88 0.035], ...
            'Callback', @(~, ~) mark_two_stage_pending());

        [state.rlamSlider, state.rlamEdit] = add_numeric_control(parent, ...
            'rlam (press Draw/Recompute)', params.rlam, 0.5, 20, 0.65);
        [state.epsSlider, state.epsEdit] = add_numeric_control(parent, ...
            'c_eps', params.c_eps, 0.05, 3, 0.56);
        [state.distSlider, state.distEdit] = add_numeric_control(parent, ...
            'c_rad', params.c_rad, 0.1, 10, 0.43);
        [state.tauSlider, state.tauEdit] = add_numeric_control(parent, ...
            'c_tau', params.c_tau, 0.05, 5, 0.30);

        state.showEdges = uicontrol(parent, 'Style', 'checkbox', ...
            'String', 'show selected edges', 'Value', 1, ...
            'Units', 'normalized', 'Position', [0.06 0.275 0.88 0.035], ...
            'Callback', @(~, ~) refresh_plot());
        state.showScaffold = uicontrol(parent, 'Style', 'checkbox', ...
            'String', 'show flat scaffold wireframe', 'Value', 0, ...
            'Units', 'normalized', 'Position', [0.06 0.24 0.88 0.035], ...
            'Callback', @(~, ~) refresh_plot());

        uicontrol(parent, 'Style', 'text', 'String', 'Color mode', ...
            'HorizontalAlignment', 'left', 'Units', 'normalized', ...
            'Position', [0.06 0.205 0.88 0.03]);
        state.colorPopup = uicontrol(parent, 'Style', 'popupmenu', ...
            'String', {'beta', 'signed mean curvature', ...
            'log10 abs mean curvature', 'patch max abs curvature', ...
            'jacobian ratio to CAD', 'displacement to CAD', 'constant'}, ...
            'Value', 3, 'Units', 'normalized', ...
            'Position', [0.06 0.175 0.88 0.04], ...
            'Callback', @(~, ~) refresh_plot());
        uicontrol(parent, 'Style', 'pushbutton', 'String', 'Update blend', ...
            'Units', 'normalized', 'Position', [0.06 0.125 0.40 0.045], ...
            'Callback', @(~, ~) apply_from_controls());
        uicontrol(parent, 'Style', 'pushbutton', ...
            'String', 'Draw/Recompute', ...
            'Units', 'normalized', 'Position', [0.54 0.125 0.40 0.045], ...
            'Callback', @(~, ~) draw_recompute());
        uicontrol(parent, 'Style', 'pushbutton', 'String', 'Save .go3', ...
            'Units', 'normalized', 'Position', [0.06 0.075 0.40 0.045], ...
            'Callback', @(~, ~) save_current_surface());
        uicontrol(parent, 'Style', 'pushbutton', 'String', 'Reset', ...
            'Units', 'normalized', 'Position', [0.54 0.075 0.40 0.045], ...
            'Callback', @(~, ~) reset_controls());
        state.status = uicontrol(parent, 'Style', 'text', ...
            'String', '', 'HorizontalAlignment', 'left', ...
            'Units', 'normalized', 'Position', [0.06 0.005 0.88 0.065]);
    end

    function [slider, editBox] = add_numeric_control(parent, label, ...
            value, minValue, maxValue, ypos)
        uicontrol(parent, 'Style', 'text', 'String', label, ...
            'HorizontalAlignment', 'left', 'Units', 'normalized', ...
            'Position', [0.06 ypos+0.08 0.88 0.04]);
        slider = uicontrol(parent, 'Style', 'slider', ...
            'Min', minValue, 'Max', maxValue, ...
            'Value', min(max(value, minValue), maxValue), ...
            'Units', 'normalized', 'Position', [0.06 ypos+0.035 0.58 0.04], ...
            'Callback', @(src, ~) slider_changed(src));
        editBox = uicontrol(parent, 'Style', 'edit', ...
            'String', sprintf('%.6g', value), ...
            'Units', 'normalized', 'Position', [0.68 ypos+0.03 0.26 0.05], ...
            'Callback', @(src, ~) edit_changed(src));
    end

    function slider_changed(src)
        if src == state.rlamSlider
            set(state.rlamEdit, 'String', sprintf('%.6g', get(src, 'Value')));
            mark_rlam_pending();
            return
        elseif src == state.epsSlider
            set(state.epsEdit, 'String', sprintf('%.6g', get(src, 'Value')));
        elseif src == state.distSlider
            set(state.distEdit, 'String', sprintf('%.6g', get(src, 'Value')));
        elseif src == state.tauSlider
            set(state.tauEdit, 'String', sprintf('%.6g', get(src, 'Value')));
        end
        apply_from_controls();
    end

    function edit_changed(src)
        if src == state.rlamEdit
            mark_rlam_pending();
        else
            apply_from_controls();
        end
    end

    function display_mode_changed()
        switch get(state.displayPopup, 'Value')
            case 1
                state.displayMode = 'blend';
            case 2
                state.displayMode = 'smooth';
            otherwise
                state.displayMode = 'scaffold';
        end
        refresh_plot();
    end

    function mark_rlam_pending()
        try
            state.pendingRlam = parse_positive_number( ...
                get(state.rlamEdit, 'String'), 'rlam');
            set(state.rlamSlider, 'Value', min(max(state.pendingRlam, ...
                get(state.rlamSlider, 'Min')), get(state.rlamSlider, 'Max')));
            set_status_text('rlam changed; press Draw/Recompute to update the surface.');
        catch ME
            set(state.status, 'String', ME.message);
        end
    end

    function mark_two_stage_pending()
        state.pendingTwoStageSmoother = logical(get( ...
            state.twoStageCheckbox, 'Value'));
        set_status_text(sprintf(['two-stage smoother changed to %s; ', ...
            'press Draw/Recompute to update the surface.'], ...
            on_off_text(state.pendingTwoStageSmoother)));
    end

    function reset_controls()
        set(state.edgeEdit, 'String', mat2str(initialParams.selectedEdgeIds));
        set(state.displayPopup, 'Value', 1);
        set(state.colorPopup, 'Value', 3);
        set(state.showEdges, 'Value', 1);
        set(state.showScaffold, 'Value', 0);
        set(state.twoStageCheckbox, 'Value', ...
            double(initialParams.useTwoStageSmoother));
        state.displayMode = 'blend';
        state.pendingTwoStageSmoother = ...
            logical(initialParams.useTwoStageSmoother);
        set(state.rlamEdit, 'String', sprintf('%.6g', initialParams.rlam));
        set(state.epsEdit, 'String', sprintf('%.6g', initialParams.c_eps));
        set(state.distEdit, 'String', sprintf('%.6g', initialParams.c_rad));
        set(state.tauEdit, 'String', sprintf('%.6g', initialParams.c_tau));
        set(state.rlamSlider, 'Value', min(max(initialParams.rlam, ...
            get(state.rlamSlider, 'Min')), get(state.rlamSlider, 'Max')));
        set(state.epsSlider, 'Value', min(max(initialParams.c_eps, ...
            get(state.epsSlider, 'Min')), get(state.epsSlider, 'Max')));
        set(state.distSlider, 'Value', min(max(initialParams.c_rad, ...
            get(state.distSlider, 'Min')), get(state.distSlider, 'Max')));
        set(state.tauSlider, 'Value', min(max(initialParams.c_tau, ...
            get(state.tauSlider, 'Min')), get(state.tauSlider, 'Max')));
        apply_from_controls();
    end

    function ok = apply_from_controls()
        ok = false;
        try
            newParams = params;
            newParams.selectedEdgeIds = parse_edge_ids( ...
                get(state.edgeEdit, 'String'), [edges.edge_id]);
            newParams.c_eps = parse_positive_number( ...
                get(state.epsEdit, 'String'), 'c_eps');
            newParams.c_rad = parse_positive_number( ...
                get(state.distEdit, 'String'), 'c_rad');
            newParams.c_tau = parse_positive_number( ...
                get(state.tauEdit, 'String'), 'c_tau');
            state.pendingRlam = parse_positive_number( ...
                get(state.rlamEdit, 'String'), 'rlam');
            state.pendingTwoStageSmoother = logical(get( ...
                state.twoStageCheckbox, 'Value'));
        catch ME
            set(state.status, 'String', ME.message);
            return
        end

        set(state.rlamSlider, 'Value', min(max(state.pendingRlam, ...
            get(state.rlamSlider, 'Min')), get(state.rlamSlider, 'Max')));
        set(state.epsSlider, 'Value', min(max(newParams.c_eps, ...
            get(state.epsSlider, 'Min')), get(state.epsSlider, 'Max')));
        set(state.distSlider, 'Value', min(max(newParams.c_rad, ...
            get(state.distSlider, 'Min')), get(state.distSlider, 'Max')));
        set(state.tauSlider, 'Value', min(max(newParams.c_tau, ...
            get(state.tauSlider, 'Min')), get(state.tauSlider, 'Max')));

        params = newParams;
        refresh_plot();
        ok = true;
    end

    function draw_recompute()
        if ~apply_from_controls()
            return
        end
        newRlam = state.pendingRlam;
        newTwoStageSmoother = state.pendingTwoStageSmoother;

        set_status_text(sprintf(['Recomputing smoother for rlam = %.6g, ', ...
            'two-stage = %s ...'], newRlam, ...
            on_off_text(newTwoStageSmoother)));
        drawnow;

        try
            optsDraw = opts;
            optsDraw.two_stage_smoother = newTwoStageSmoother;
            [state.S_smooth, sigmaCadNew, gradSigmaCadNew] = ...
                recompute_smooth_state(scaffoldFile, norder, optsDraw, ...
                S_cad, newRlam);
            state.reportSmooth = surface_validity_report(state.S_smooth, ...
                'S_smooth_live', S_cad);
            print_surface_validity_report(state.reportSmooth);
            params.rlam = newRlam;
            params.useTwoStageSmoother = newTwoStageSmoother;
            params.sigma = sigmaCadNew;
            params.gradSigma = gradSigmaCadNew;
            state.S_blend = [];
            state.reportBlend = [];
            state.axisLimits = geometry_axis_limits(S_cad, ...
                state.S_smooth, scaffold);
            refresh_plot();
        catch ME
            set(state.status, 'String', ME.message);
        end
    end

    function refresh_plot()
        if strcmp(state.displayMode, 'blend')
            try
                [state.S_blend, state.beta, state.dsoft] = ...
                    build_blended_surface_spectral(S_cad, state.S_smooth, ...
                    edges, params);
                state.reportBlend = surface_validity_report(state.S_blend, ...
                    'S_blend_live', S_cad);
            catch ME
                set(state.status, 'String', ME.message);
                return
            end
        end

        camera = save_camera(ax);
        delete_graphics(state.surfaceHandle);
        delete_graphics(state.edgeHandles);
        delete_graphics(state.scaffoldHandle);
        cla(ax);

        axes(ax);
        if strcmp(state.displayMode, 'scaffold')
            state.surfaceHandle = plot_flat_scaffold(ax, scaffold, 'solid');
            state.edgeHandles = gobjects(0);
            state.scaffoldHandle = gobjects(0);
            clim(ax, [0 1]);
            cb = colorbar(ax);
            ylabel(cb, 'flat scaffold', 'Interpreter', 'none');
            apply_fixed_axis_limits(ax, state.axisLimits);
            view(ax, 3);
            camlight headlight;
            lighting(ax, 'gouraud');
            grid(ax, 'on');
            xlabel(ax, 'x');
            ylabel(ax, 'y');
            zlabel(ax, 'z');
            title(ax, [packageName ': flat .gidmsh scaffold'], ...
                'Interpreter', 'none');
            restore_camera(ax, camera);
            set_status_text(sprintf(['mode flat scaffold\n', ...
                'triangles %d, vertices %d\n', ...
                'flat scaffold: no high-order curvature'], ...
                size(scaffold.tris, 2), size(scaffold.nodes, 2)));
            drawnow limitrate;
            return
        end

        colorMode = selected_color_mode();
        if strcmp(state.displayMode, 'blend')
            S_display = state.S_blend;
            reportDisplay = state.reportBlend;
            displayTitle = 'live edge-preserve blend';
        else
            S_display = state.S_smooth;
            reportDisplay = state.reportSmooth;
            displayTitle = 'fully smooth CAD-skeleton surface';
        end

        betaForColor = state.beta;
        if strcmp(state.displayMode, 'smooth')
            betaForColor = [];
        end
        [plotValues, colorLabel, colorStatusNote] = display_color_values( ...
            S_display, S_cad, betaForColor, colorMode);
        if strcmp(colorMode, 'constant') || ...
                strcmp(colorStatusNote, 'beta unavailable in fully smooth mode')
            plotArgsUse = [surfacePlotArgs, ...
                {'FaceColor', [0.55 0.75 0.95]}];
        else
            plotArgsUse = surfacePlotArgs;
        end
        state.surfaceHandle = plot_surfer_smooth_lit(S_display, ...
            plotValues, plotArgsUse, surfacePlotRefine);

        hold(ax, 'on');
        if logical(get(state.showEdges, 'Value'))
            state.edgeHandles = plot_selected_edges(ax, edges, ...
                params.selectedEdgeIds, [0 0 0], 3.0);
        else
            state.edgeHandles = gobjects(0);
        end
        if logical(get(state.showScaffold, 'Value'))
            state.scaffoldHandle = plot_flat_scaffold(ax, scaffold, ...
                'wireframe');
        else
            state.scaffoldHandle = gobjects(0);
        end
        apply_fixed_axis_limits(ax, state.axisLimits);
        view(ax, 3);
        camlight headlight;
        colorLimitText = apply_color_limits(ax, plotValues, colorMode, ...
            colorStatusNote);
        cb = colorbar(ax);
        ylabel(cb, colorLabel, 'Interpreter', 'none');
        grid(ax, 'on');
        xlabel(ax, 'x');
        ylabel(ax, 'y');
        zlabel(ax, 'z');
        title(ax, sprintf('%s: %s | color: %s', packageName, ...
            displayTitle, colorLabel), 'Interpreter', 'none');
        restore_camera(ax, camera);

        healthText = report_health_text(reportDisplay);
        if strcmp(state.displayMode, 'blend')
            set_status_text(sprintf([ ...
                'mode blend | edges %s\n', ...
                'color %s%s%s\n', ...
                'rlam drawn/pending %.4g %.4g\n', ...
                'two-stage drawn/pending %s %s\n', ...
                'c eps/rad/tau %.4g %.4g %.4g\n', ...
                'sigma min/mean/max %.3g %.3g %.3g\n', ...
                'd min/mean/max %.3g %.3g %.3g\n', ...
                'b(0 CAD,1 smooth) min/mean/max %.3g %.3g %.3g\n', ...
                '%s'], ...
                mat2str(params.selectedEdgeIds), colorLabel, ...
                status_note_suffix(colorStatusNote), colorLimitText, ...
                params.rlam, state.pendingRlam, ...
                on_off_text(params.useTwoStageSmoother), ...
                on_off_text(state.pendingTwoStageSmoother), ...
                params.c_eps, params.c_rad, params.c_tau, ...
                min(local_sigma(params)), ...
                mean(local_sigma(params)), max(local_sigma(params)), ...
                min(state.dsoft), mean(state.dsoft), max(state.dsoft), ...
                min(state.beta), mean(state.beta), max(state.beta), ...
                healthText));
        else
            set_status_text(sprintf([ ...
                'mode fully smooth\n', ...
                'color %s%s%s\n', ...
                'rlam drawn/pending %.4g %.4g\n', ...
                'two-stage drawn/pending %s %s\n', ...
                'shown surface ignores edge/beta controls\n', ...
                'sigma min/mean/max %.3g %.3g %.3g\n', ...
                '%s'], ...
                colorLabel, status_note_suffix(colorStatusNote), ...
                colorLimitText, params.rlam, state.pendingRlam, ...
                on_off_text(params.useTwoStageSmoother), ...
                on_off_text(state.pendingTwoStageSmoother), ...
                min(local_sigma(params)), mean(local_sigma(params)), ...
                max(local_sigma(params)), healthText));
        end
        drawnow limitrate;
    end

    function save_current_surface()
        if ~apply_from_controls()
            return
        end

        [S_to_save, saveMode] = current_display_surface();
        if isempty(S_to_save)
            set_status_text('No displayed surface is available to save yet.');
            return
        end

        outputDir = fullfile(repoRoot, 'geometries', 'meshes', ...
            'edge_preserve_outputs', packageName);
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end
        suggestedName = suggested_surface_filename(packageName, params, saveMode);
        [fname, fpath] = uiputfile({'*.go3', 'GO3 surface (*.go3)'}, ...
            'Save displayed surface as .go3', fullfile(outputDir, suggestedName));
        if isequal(fname, 0) || isequal(fpath, 0)
            set_status_text('Save cancelled.');
            return
        end
        if ~endsWith(fname, '.go3', 'IgnoreCase', true)
            fname = [fname '.go3'];
        end
        outFile = fullfile(fpath, fname);

        try
            reportSave = surface_validity_report(S_to_save, ...
                ['save_' saveMode], S_cad);
            print_surface_validity_report(reportSave);
            assert_surface_valid_for_bie(reportSave, allowInvalidGo3Save);

            write_go3_file(outFile, S_to_save);
            S_saved = surfsmooth3d.surfer.load_from_file(outFile);
            assert_compatible_surfaces(S_cad, S_saved);
            reportSaved = surface_validity_report(S_saved, ...
                ['saved_' saveMode], S_cad);
            print_surface_validity_report(reportSaved);
            assert_surface_valid_for_bie(reportSaved, allowInvalidGo3Save);
            normalDotMin = min(sum(S_saved.n .* S_cad.n, 1));
            pendingText = '';
            if abs(state.pendingRlam - params.rlam) > ...
                    10*eps(max(1, abs(params.rlam)))
                pendingText = sprintf(['\nNote: pending rlam %.6g was ', ...
                    'not drawn; saved drawn rlam %.6g.'], ...
                    state.pendingRlam, params.rlam);
            end
            if state.pendingTwoStageSmoother ~= ...
                    logical(params.useTwoStageSmoother)
                pendingText = [pendingText, sprintf([ ...
                    '\nNote: pending two-stage smoother %s was not drawn; ', ...
                    'saved drawn two-stage smoother %s.'], ...
                    on_off_text(state.pendingTwoStageSmoother), ...
                    on_off_text(params.useTwoStageSmoother))];
            end
            set_status_text(sprintf(['Saved %s .go3:\n%s\n', ...
                'min saved-normal dot CAD-normal %.3g%s'], ...
                saveMode, outFile, normalDotMin, pendingText));
            fprintf('Saved %s .go3: %s\n', saveMode, outFile);
        catch ME
            set(state.status, 'String', ME.message);
        end
    end

    function [S_display, saveMode] = current_display_surface()
        if strcmp(state.displayMode, 'blend')
            if isempty(state.S_blend)
                try
                    [state.S_blend, state.beta, state.dsoft] = ...
                        build_blended_surface_spectral(S_cad, ...
                        state.S_smooth, edges, params);
                catch ME
                    set(state.status, 'String', ME.message);
                    S_display = [];
                    saveMode = 'edgepreserve';
                    return
                end
            end
            S_display = state.S_blend;
            saveMode = 'edgepreserve';
        else
            S_display = state.S_smooth;
            saveMode = 'fullysmooth';
        end
        if strcmp(state.displayMode, 'scaffold')
            set_status_text('The flat scaffold is a visualization only; save a smooth or blend .go3 surface instead.');
            S_display = [];
            saveMode = 'scaffold';
        end
    end

    function colorMode = selected_color_mode()
        entries = get(state.colorPopup, 'String');
        colorMode = entries{get(state.colorPopup, 'Value')};
    end

    function set_status_text(txt)
        set(state.status, 'String', txt);
    end
end

function [plotValues, colorLabel, statusNote] = display_color_values( ...
        S_display, S_cad, beta, colorMode)
    statusNote = '';
    switch colorMode
        case 'beta'
            if numel(beta) == S_display.npts
                plotValues = beta(:);
                colorLabel = 'beta (0 CAD, 1 smooth)';
            else
                plotValues = zeros(S_display.npts, 1);
                colorLabel = 'constant';
                statusNote = 'beta unavailable in fully smooth mode';
            end
        case 'signed mean curvature'
            plotValues = S_display.mean_curv(:);
            colorLabel = 'signed mean curvature';
        case 'log10 abs mean curvature'
            plotValues = log10(abs(S_display.mean_curv(:)) + 1e-16);
            colorLabel = 'log10(|mean curvature| + 1e-16)';
        case 'patch max abs curvature'
            plotValues = patch_max_abs_curvature_values(S_display);
            colorLabel = 'patch max |mean curvature|';
        case 'jacobian ratio to CAD'
            [srcDisplay, ~, ~, ~, ~, ~] = extract_arrays(S_display);
            [srcCad, ~, ~, ~, ~, ~] = extract_arrays(S_cad);
            jacDisplay = surface_jacobian_from_srcvals(srcDisplay);
            jacCad = surface_jacobian_from_srcvals(srcCad);
            plotValues = jacDisplay ./ max(jacCad, realmin);
            colorLabel = '|du x dv| / CAD';
        case 'displacement to CAD'
            plotValues = sqrt(sum((S_display.r - S_cad.r).^2, 1)).';
            colorLabel = 'displacement to CAD';
        case 'constant'
            plotValues = zeros(S_display.npts, 1);
            colorLabel = 'constant';
        otherwise
            error('driver_edge_preserve_compare_surfaces:badColorMode', ...
                'Unknown color mode "%s".', colorMode);
    end
end

function limitText = apply_color_limits(ax, plotValues, colorMode, statusNote)
    if strcmp(colorMode, 'signed mean curvature')
        finiteVals = plotValues(isfinite(plotValues));
        if isempty(finiteVals)
            climVals = [-1 1];
            statsText = 'signed H unavailable';
        else
            maxAbs = max(abs(finiteVals));
            if maxAbs <= 0
                maxAbs = 1;
            end
            climVals = [-maxAbs maxAbs];
            statsText = sprintf('signed H min/mean/max %.3g %.3g %.3g', ...
                min(finiteVals), mean(finiteVals), max(finiteVals));
        end
        clim(ax, climVals);
        limitText = sprintf('\n%s | clim [%.3g %.3g]', ...
            statsText, climVals(1), climVals(2));
        return
    end

    if strcmp(colorMode, 'beta') && isempty(statusNote)
        clim(ax, [0 1]);
        limitText = sprintf('\nclim [0 1]');
        return
    end

    if strcmp(colorMode, 'constant') || ...
            strcmp(statusNote, 'beta unavailable in fully smooth mode')
        clim(ax, [0 1]);
        limitText = sprintf('\nclim [0 1]');
        return
    end

    finiteVals = plotValues(isfinite(plotValues));
    if isempty(finiteVals)
        clim(ax, [0 1]);
        limitText = sprintf('\nclim [0 1]');
        return
    end

    lo = min(finiteVals);
    hi = max(finiteVals);
    if lo == hi
        pad = max(1, abs(lo))*1e-6;
        lo = lo - pad;
        hi = hi + pad;
    end
    clim(ax, [lo hi]);
    limitText = sprintf('\nclim [%.3g %.3g]', lo, hi);
end

function values = patch_max_abs_curvature_values(S)
    values = accumarray(S.patch_id(:), abs(S.mean_curv(:)), ...
        [S.npatches 1], @max, nan);
end

function jac = surface_jacobian_from_srcvals(srcvals)
    cr = cross(srcvals(4:6, :).', srcvals(7:9, :).').';
    jac = sqrt(sum(cr.^2, 1)).';
end

function txt = report_health_text(report)
    if isempty(report)
        txt = 'validity unavailable';
        return
    end
    if ~report.is_valid
        status = 'FAIL';
    elseif report.has_warnings
        status = 'WARN';
    else
        status = 'PASS';
    end

    [worstPatch, worstPatchCurv] = worst_patch_by_abs_curvature(report);
    txt = sprintf(['validity %s | max |H| %.3g | worst patch %d %.3g\n', ...
        'min jac/ref %.3g | bad nodes/patches %d/%d'], ...
        status, report.stats.mean_curv_abs_max, worstPatch, ...
        worstPatchCurv, report.stats.jac_ratio_min, ...
        report.stats.bad_node_count, report.stats.bad_patch_count);
end

function [patchIndex, patchCurv] = worst_patch_by_abs_curvature(report)
    patchCurvVals = accumarray(report.patch_id(:), ...
        abs(report.mean_curv(:)), [report.stats.npatches 1], @max, nan);
    [patchCurv, patchIndex] = max(patchCurvVals);
end

function suffix = status_note_suffix(statusNote)
    if isempty(statusNote)
        suffix = '';
    else
        suffix = sprintf(' (%s)', statusNote);
    end
end

function limits = geometry_axis_limits(S_cad, S_smooth, scaffold)
    xyz = [S_cad.r, S_smooth.r, scaffold.nodes];
    mins = min(xyz, [], 2);
    maxs = max(xyz, [], 2);
    center = 0.5 * (mins + maxs);
    span = max(maxs - mins);
    span = max(span, 1);
    halfWidth = 0.55 * span;
    limits = [center(1)-halfWidth, center(1)+halfWidth; ...
              center(2)-halfWidth, center(2)+halfWidth; ...
              center(3)-halfWidth, center(3)+halfWidth];
end

function apply_fixed_axis_limits(ax, limits)
    axis(ax, 'equal');
    xlim(ax, limits(1, :));
    ylim(ax, limits(2, :));
    zlim(ax, limits(3, :));
    axis(ax, 'vis3d');
end

function camera = save_camera(ax)
    camera = struct();
    if isempty(ax.Children)
        camera.valid = false;
        return
    end
    camera.valid = true;
    camera.pos = campos(ax);
    camera.target = camtarget(ax);
    camera.up = camup(ax);
    camera.va = camva(ax);
end

function restore_camera(ax, camera)
    if ~camera.valid
        return
    end
    campos(ax, camera.pos);
    camtarget(ax, camera.target);
    camup(ax, camera.up);
    camva(ax, camera.va);
end

function delete_graphics(h)
    h = h(isgraphics(h));
    if ~isempty(h)
        delete(h);
    end
end

function write_go3_file(outFile, S)
    [srcvals, ~, ~, ~, ~, ~] = extract_arrays(S);
    norder = infer_common_order(S);

    fid = fopen(outFile, 'w');
    if fid < 0
        error('driver_edge_preserve_compare_surfaces:go3OpenFailed', ...
            'Could not open output .go3 file for writing: %s', outFile);
    end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '%d\n', norder);
    fprintf(fid, '%d\n', S.npatches);
    fprintf(fid, '%.16e\n', srcvals.');
    clear cleanup;
end

function fname = suggested_surface_filename(packageName, params, saveMode)
    smootherText = smoother_stage_filename_text(params);
    switch saveMode
        case 'edgepreserve'
            edgeText = strjoin(arrayfun(@(id) sprintf('%d', id), ...
                params.selectedEdgeIds(:).', 'UniformOutput', false), '-');
            if isempty(edgeText)
                edgeText = 'none';
            end
            fname = sprintf(['%s_edgepreserve_%s_edges%s_rlam%s_', ...
                'ceps%s_crad%s_ctau%s.go3'], ...
                sanitize_filename_text(packageName), smootherText, ...
                sanitize_filename_text(edgeText), ...
                sanitize_number_for_filename(params.rlam), ...
                sanitize_number_for_filename(params.c_eps), ...
                sanitize_number_for_filename(params.c_rad), ...
                sanitize_number_for_filename(params.c_tau));
        case 'fullysmooth'
            fname = sprintf('%s_fullysmooth_%s_rlam%s.go3', ...
                sanitize_filename_text(packageName), smootherText, ...
                sanitize_number_for_filename(params.rlam));
        otherwise
            error('driver_edge_preserve_compare_surfaces:badSaveMode', ...
                'Unknown save mode "%s".', saveMode);
    end
end

function txt = smoother_stage_filename_text(params)
    if isfield(params, 'useTwoStageSmoother') && ...
            logical(params.useTwoStageSmoother)
        txt = 'twostage';
    else
        txt = 'onestage';
    end
end

function txt = on_off_text(tf)
    if logical(tf)
        txt = 'on';
    else
        txt = 'off';
    end
end

function txt = sanitize_number_for_filename(x)
    txt = sprintf('%.6g', x);
    txt = strrep(txt, '-', 'm');
    txt = strrep(txt, '+', '');
    txt = strrep(txt, '.', 'p');
    txt = strrep(txt, 'e', 'e');
    txt = sanitize_filename_text(txt);
end

function txt = sanitize_filename_text(txt)
    txt = char(txt);
    txt = regexprep(txt, '[^A-Za-z0-9_-]+', '_');
    txt = regexprep(txt, '_+', '_');
    txt = regexprep(txt, '^_+|_+$', '');
    if isempty(txt)
        txt = 'surface';
    end
end

function val = parse_positive_number(str, name)
    val = str2double(str);
    if ~isfinite(val) || val <= 0
        error('driver_edge_preserve_compare_surfaces:badParameter', ...
            '%s must be a positive finite scalar.', name);
    end
end

function ids = parse_edge_ids(str, validIds)
    tokens = regexp(str, '[-+]?\d+', 'match');
    ids = unique(str2double(tokens));
    ids = ids(isfinite(ids));
    ids = ids(ismember(ids, validIds));
    if isempty(ids)
        error('driver_edge_preserve_compare_surfaces:noValidEdgeIds', ...
            'No valid CAD edge IDs were selected.');
    end
end

function edges = read_edges_v2(fid, edgeFile)
    edgeOrder = read_scalar_line(fid, 'edge_order');
    numEdges = read_scalar_line(fid, 'num_cad_edges');
    edges = repmat(empty_edge(edgeOrder), 1, numEdges);
    for k = 1:numEdges
        expect_line(fid, 'CAD_EDGE', edgeFile);
        edges(k).edge_order = edgeOrder;
        edges(k).format_version = 2;
        edges(k).edge_id = read_scalar_line(fid, 'edge_id');
        edges(k).curve_tag = read_scalar_line(fid, 'curve_tag');
        edges(k).adjacent_faces = read_vector_line(fid, 'adjacent_faces');
        edges(k).num_panels = read_scalar_line(fid, 'num_panels');
        edges(k).panels = repmat(empty_panel(edgeOrder), ...
            1, edges(k).num_panels);
        for j = 1:edges(k).num_panels
            edges(k).panels(j) = read_edge_panel(fid, edgeFile, edgeOrder);
        end
        expect_line(fid, 'END_CAD_EDGE', edgeFile);
    end
    expect_line(fid, 'END_STEP_MESHER_EDGES_GLL', edgeFile);
end

function edges = read_edges_v1_as_v2(fid, edgeFile)
    edgeOrder = read_scalar_line(fid, 'edge_order');
    numEdges = read_scalar_line(fid, 'num_edges');
    edges = repmat(empty_edge(edgeOrder), 1, numEdges);
    for k = 1:numEdges
        expect_line(fid, 'EDGE', edgeFile);
        edges(k).edge_order = edgeOrder;
        edges(k).format_version = 1;
        edges(k).edge_id = read_scalar_line(fid, 'edge_id');
        edges(k).curve_tag = read_scalar_line(fid, 'curve_tag');
        nodeTags = read_vector_line(fid, 'node_tags');
        edges(k).adjacent_faces = read_vector_line(fid, 'adjacent_faces');
        numNodes = read_scalar_line(fid, 'num_nodes');
        line = next_nonempty_line(fid);
        if ~startsWith(line, 'columns ')
            error('driver_edge_preserve_compare_surfaces:badEdgesFile', ...
                'Expected columns line in %s.', edgeFile);
        end
        data = fscanf(fid, '%f', [8, numNodes]).';
        panel = empty_panel(edgeOrder);
        panel.panel_id = 1;
        panel.node_tags = nodeTags;
        panel.num_nodes = numNodes;
        panel.a = data(:, 1).';
        panel.xyz = data(:, 2:4).';
        panel.tangent = data(:, 5:7).';
        panel.w = data(:, 8).';
        edges(k).num_panels = 1;
        edges(k).panels = panel;
        expect_line(fid, 'END_EDGE', edgeFile);
    end
    expect_line(fid, 'END_STEP_MESHER_EDGES_GLL', edgeFile);
end

function panel = read_edge_panel(fid, edgeFile, edgeOrder)
    expect_line(fid, 'PANEL', edgeFile);
    panel = empty_panel(edgeOrder);
    panel.panel_id = read_scalar_line(fid, 'panel_id');
    panel.node_tags = read_vector_line(fid, 'node_tags');
    panel.num_nodes = read_scalar_line(fid, 'num_nodes');
    line = next_nonempty_line(fid);
    if ~startsWith(line, 'columns ')
        error('driver_edge_preserve_compare_surfaces:badEdgesFile', ...
            'Expected columns line in %s.', edgeFile);
    end
    data = fscanf(fid, '%f', [8, panel.num_nodes]).';
    panel.a = data(:, 1).';
    panel.xyz = data(:, 2:4).';
    panel.tangent = data(:, 5:7).';
    panel.w = data(:, 8).';
    expect_line(fid, 'END_PANEL', edgeFile);
end

function edge = empty_edge(edgeOrder)
    edge = struct();
    edge.edge_order = edgeOrder;
    edge.format_version = [];
    edge.edge_id = [];
    edge.curve_tag = [];
    edge.adjacent_faces = [];
    edge.num_panels = [];
    edge.panels = repmat(empty_panel(edgeOrder), 1, 0);
end

function panel = empty_panel(edgeOrder)
    panel = struct();
    panel.edge_order = edgeOrder;
    panel.panel_id = [];
    panel.node_tags = [];
    panel.num_nodes = [];
    panel.a = [];
    panel.xyz = [];
    panel.tangent = [];
    panel.w = [];
end

function val = read_scalar_line(fid, key)
    vals = read_vector_line(fid, key);
    val = vals(1);
end

function vals = read_vector_line(fid, key)
    line = next_nonempty_line(fid);
    parts = split(strtrim(line));
    if isempty(parts) || ~strcmp(parts{1}, key)
        error('driver_edge_preserve_compare_surfaces:badEdgesFile', ...
            'Expected line beginning with "%s".', key);
    end
    vals = str2double(parts(2:end)).';
end

function expect_line(fid, expected, edgeFile)
    line = next_nonempty_line(fid);
    if ~strcmp(line, expected)
        error('driver_edge_preserve_compare_surfaces:badEdgesFile', ...
            'Expected "%s" in %s, got "%s".', expected, edgeFile, line);
    end
end

function line = next_nonempty_line(fid)
    while true
        line = fgetl(fid);
        if ~ischar(line)
            error('driver_edge_preserve_compare_surfaces:unexpectedEof', ...
                'Unexpected end of edge file.');
        end
        line = strtrim(line);
        if ~isempty(line)
            return
        end
    end
end

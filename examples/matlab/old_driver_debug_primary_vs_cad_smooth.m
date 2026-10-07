%DRIVER_DEBUG_PRIMARY_VS_CAD_SMOOTH Compare primary and CAD-skeleton smoothers.
%
% This diagnostic runs the original surface smoother with no CAD skeleton
% on the same scaffold.gidmsh used by the CAD-skeleton smoother, then
% compares it against the CAD-point-skeleton smoothed surface. The goal is
% to decide whether bad .go3 patches are caused by the base scaffold /
% smoother machinery or specifically by the CAD-skeleton Phi path.

clear;
format long g;

driverFile = mfilename('fullpath');
examplesDir = fileparts(driverFile);
repoRoot = fileparts(fileparts(examplesDir));
matlabDir = fullfile(repoRoot, 'matlab');
addpath(matlabDir);
setup_surfsmooth3d;
packageRoot = surfsmooth3d.edgepreserve.default_input_root(repoRoot);

%% User parameters
packageDir = fullfile(repoRoot, 'tests', 'fixtures', 'cad', ...
    'convergence_intersection_two_balls', 'p4_mf0p10000_rl0');

% Optional existing CAD-skeleton smoother output. The default recomputes
% it from surface.go3. Set this path and cadSmoothSource = 'saved' to compare
% against a previously exported surface with the same patch layout.
savedCadSmoothFile = '';

rlam = 2;
patchIndex = 1;
cadSmoothSource = default_cad_smooth_source();
% cadSmoothSource = 'saved';
nWorst = 20;
denseRefine = 34;
makePlots = usejava('desktop');

%% Load package inputs
cadFile = fullfile(packageDir, 'surface.go3');
scaffoldFile = fullfile(packageDir, 'scaffold.gidmsh');

if ~exist(cadFile, 'file')
    error('driver_debug_primary_vs_cad_smooth:missingCadFile', ...
        'CAD surface file not found:\n%s', cadFile);
end
if ~exist(scaffoldFile, 'file')
    error('driver_debug_primary_vs_cad_smooth:missingScaffoldFile', ...
        'Scaffold mesh file not found:\n%s', scaffoldFile);
end

S_cad = surfsmooth3d.surfer.load_from_file(cadFile);
norder = infer_common_order(S_cad);

fprintf('driver_debug_primary_vs_cad_smooth\n');
fprintf('  packageDir: %s\n', packageDir);
fprintf('  scaffold: %s\n', scaffoldFile);
fprintf('  CAD surface: %s\n', cadFile);
fprintf('  norder=%d, rlam=%g, patchIndex=%d\n', ...
    norder, rlam, patchIndex);

%% Run original/no-CAD primary smoother on the same scaffold
optsPrimary = struct();
optsPrimary.filetype = 3;
optsPrimary.nrefine = 0;
optsPrimary.nquad = norder;
optsPrimary.rlam = rlam;

fprintf('\nRunning primary smoother: no opts.fcad, same scaffold...\n');
tic;
S_primary_all = surfsmooth3d.multiscale_mesher(scaffoldFile, norder, optsPrimary);
primaryTime = toc;
S_primary = S_primary_all{end};
fprintf('  primary smoother time = %.3f s\n', primaryTime);

%% Load or recompute CAD-point-skeleton smoother result
if strcmp(cadSmoothSource, 'recompute')
    fprintf('\nRunning CAD-point-skeleton smoother from surface.go3...\n');
    tempDir = fullfile(tempdir, 'surfsmooth3d_primary_vs_cad_smooth');
    if ~exist(tempDir, 'dir')
        mkdir(tempDir);
    end
    [~, packageName] = fileparts(packageDir);
    tempCadSkeletonTxt = fullfile(tempDir, ...
        sprintf('%s_cad_skeleton_from_go3.txt', packageName));
    write_cad_skeleton_txt(S_cad, tempCadSkeletonTxt);

    optsCad = optsPrimary;
    optsCad.fcad = tempCadSkeletonTxt;
    tic;
    S_cad_smooth_all = surfsmooth3d.multiscale_mesher(scaffoldFile, norder, optsCad);
    cadSmoothTime = toc;
    S_cad_smooth = S_cad_smooth_all{end};
    fprintf('  CAD-point smoother time = %.3f s\n', cadSmoothTime);
    fprintf('  CAD skeleton txt: %s\n', tempCadSkeletonTxt);
elseif strcmp(cadSmoothSource, 'saved')
    if ~exist(savedCadSmoothFile, 'file')
        error('driver_debug_primary_vs_cad_smooth:missingSavedCadSmooth', ...
            ['Saved CAD-smooth file not found:\n%s\n', ...
            'Set cadSmoothSource = ''recompute'' to regenerate it.'], ...
            savedCadSmoothFile);
    end
    fprintf('\nLoading saved CAD-point-skeleton smoother output:\n  %s\n', ...
        savedCadSmoothFile);
    S_cad_smooth = surfsmooth3d.surfer.load_from_file(savedCadSmoothFile);
else
    error('driver_debug_primary_vs_cad_smooth:badCadSmoothSource', ...
        'Unknown cadSmoothSource "%s". Use "saved" or "recompute".', ...
        cadSmoothSource);
end

assert_compatible_surfaces(S_cad, S_primary, 'S_cad', 'S_primary');
assert_compatible_surfaces(S_cad, S_cad_smooth, 'S_cad', 'S_cad_smooth');

%% Global validity comparison
reportCad = surface_report(S_cad, S_cad, 'CAD input');
reportPrimary = surface_report(S_primary, S_cad, 'primary smoother');
reportCadSmooth = surface_report(S_cad_smooth, S_cad, ...
    'CAD-point skeleton smoother');

print_report(reportCad);
print_report(reportPrimary);
print_report(reportCadSmooth);

fprintf('\nWorst patches by max |mean curvature|, primary smoother:\n');
disp(worst_patch_table(S_primary, S_cad, nWorst));

fprintf('\nWorst patches by max |mean curvature|, CAD-point skeleton smoother:\n');
disp(worst_patch_table(S_cad_smooth, S_cad, nWorst));

%% Focused patch report
fprintf('\nPatch %d comparison:\n', patchIndex);
print_patch_summary(S_cad, S_cad, patchIndex, 'CAD input');
print_patch_summary(S_primary, S_cad, patchIndex, 'primary smoother');
print_patch_summary(S_cad_smooth, S_cad, patchIndex, ...
    'CAD-point skeleton smoother');

%% Plots
if makePlots
    plot_global_scalar(S_primary, log10(abs(S_primary.mean_curv(:)) + 1e-16), ...
        'primary smoother: log10(|mean curvature|)');
    plot_global_scalar(S_cad_smooth, ...
        log10(abs(S_cad_smooth.mean_curv(:)) + 1e-16), ...
        'CAD-point skeleton smoother: log10(|mean curvature|)');

    plot_global_scalar(S_primary, patch_max_abs_curvature(S_primary), ...
        'primary smoother: patch max |mean curvature|');
    plot_global_scalar(S_cad_smooth, ...
        patch_max_abs_curvature(S_cad_smooth), ...
        'CAD-point skeleton smoother: patch max |mean curvature|');

    plot_patch_overlay(S_cad, S_primary, S_cad_smooth, patchIndex, ...
        denseRefine);
    plot_patch_diagnostics(S_primary, S_cad, patchIndex, denseRefine, ...
        'primary smoother');
    plot_patch_diagnostics(S_cad_smooth, S_cad, patchIndex, denseRefine, ...
        'CAD-point skeleton smoother');
else
    fprintf('  plotting disabled because MATLAB desktop is unavailable\n');
end

function norder = infer_common_order(S)
    if isempty(S.norders)
        error('driver_debug_primary_vs_cad_smooth:emptyOrders', ...
            'Surface has no patch orders.');
    end
    if any(S.norders ~= S.norders(1))
        error('driver_debug_primary_vs_cad_smooth:mixedOrders', ...
            'This driver expects one common order across all patches.');
    end
    norder = S.norders(1);
end

function source = default_cad_smooth_source()
    source = 'recompute';
end

function write_cad_skeleton_txt(S, outFile)
    [srcvals, ~, ~, ~, ~, wts] = extract_arrays(S);
    data = [srcvals(1:3, :).', srcvals(10:12, :).', wts(:)];

    fid = fopen(outFile, 'w');
    if fid < 0
        error('driver_debug_primary_vs_cad_smooth:openFailed', ...
            'Could not open temporary CAD skeleton file: %s', outFile);
    end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '# x y z nx ny nz w generated from surface.go3\n');
    fprintf(fid, '%.16e %.16e %.16e %.16e %.16e %.16e %.16e\n', data.');
    clear cleanup;
end

function assert_compatible_surfaces(S_a, S_b, labelA, labelB)
    if S_a.npatches ~= S_b.npatches || S_a.npts ~= S_b.npts || ...
            ~isequal(S_a.norders(:), S_b.norders(:)) || ...
            ~isequal(S_a.iptype(:), S_b.iptype(:))
        error('driver_debug_primary_vs_cad_smooth:surfaceMismatch', ...
            '%s and %s do not have matching patch/node layout.', ...
            labelA, labelB);
    end
end

function report = surface_report(S, S_ref, label)
    [srcvals, ~, ~, ~, ~, wts] = extract_arrays(S);
    jac = surface_jacobian_from_srcvals(srcvals);
    [srcRef, ~, ~, ~, ~, ~] = extract_arrays(S_ref);
    jacRef = surface_jacobian_from_srcvals(srcRef);
    cr = cross(srcvals(4:6, :).', srcvals(7:9, :).').';
    crUnit = cr ./ max(jac.', realmin);
    normalCrossDot = sum(srcvals(10:12, :) .* crUnit, 1).';

    report = struct();
    report.label = label;
    report.npatches = S.npatches;
    report.npts = S.npts;
    report.area = S.area;
    report.normal_cross_dot = normalCrossDot;
    report.jac = jac;
    report.jac_ratio = jac ./ max(jacRef, realmin);
    report.disp_to_ref = sqrt(sum((S.r - S_ref.r).^2, 1)).';
    report.wts = wts(:);
    report.mean_curv = S.mean_curv(:);
    report.bad_normal = normalCrossDot < 0.99;
    report.bad_jac_ratio = report.jac_ratio < 1e-3;
    report.bad_curv = abs(report.mean_curv) > ...
        100 * max(1, max(abs(S_ref.mean_curv(:))));
    report.bad_mask = report.bad_normal | report.bad_jac_ratio | ...
        report.bad_curv | ~isfinite(report.mean_curv) | ...
        ~isfinite(report.jac) | ~isfinite(report.normal_cross_dot);
end

function print_report(report)
    fprintf('\nSurface report: %s\n', report.label);
    fprintf('  npatches=%d, npts=%d, area=%.16e\n', ...
        report.npatches, report.npts, report.area);
    fprintf('  n dot normalize(du x dv) min/mean/max = %.6e %.6e %.6e\n', ...
        min(report.normal_cross_dot), mean(report.normal_cross_dot), ...
        max(report.normal_cross_dot));
    fprintf('  |du x dv| min/mean/max = %.6e %.6e %.6e\n', ...
        min(report.jac), mean(report.jac), max(report.jac));
    fprintf('  jac/ref min/mean/max = %.6e %.6e %.6e\n', ...
        min(report.jac_ratio), mean(report.jac_ratio), ...
        max(report.jac_ratio));
    fprintf('  displacement to CAD min/mean/max/std = %.6e %.6e %.6e %.6e\n', ...
        min(report.disp_to_ref), mean(report.disp_to_ref), ...
        max(report.disp_to_ref), std(report.disp_to_ref));
    fprintf('  mean curvature min/mean/max/std = %.6e %.6e %.6e %.6e\n', ...
        min(report.mean_curv), mean(report.mean_curv), ...
        max(report.mean_curv), std(report.mean_curv));
    fprintf('  bad nodes: total=%d, normal=%d, jac_ratio=%d, curvature=%d\n', ...
        nnz(report.bad_mask), nnz(report.bad_normal), ...
        nnz(report.bad_jac_ratio), nnz(report.bad_curv));
end

function T = worst_patch_table(S, S_ref, nWorst)
    [srcvals, ~, ~, ~, ~, ~] = extract_arrays(S);
    [srcRef, ~, ~, ~, ~, ~] = extract_arrays(S_ref);
    jac = surface_jacobian_from_srcvals(srcvals);
    jacRef = surface_jacobian_from_srcvals(srcRef);
    jacRatio = jac ./ max(jacRef, realmin);
    dispToRef = sqrt(sum((S.r - S_ref.r).^2, 1)).';
    absCurv = abs(S.mean_curv(:));

    patchMaxAbsCurv = accumarray(S.patch_id(:), absCurv, ...
        [S.npatches 1], @max, nan);
    patchMeanAbsCurv = accumarray(S.patch_id(:), absCurv, ...
        [S.npatches 1], @mean, nan);
    patchMinJacRatio = accumarray(S.patch_id(:), jacRatio, ...
        [S.npatches 1], @min, nan);
    patchMaxDisp = accumarray(S.patch_id(:), dispToRef, ...
        [S.npatches 1], @max, nan);
    [~, order] = sort(patchMaxAbsCurv, 'descend');
    order = order(1:min(nWorst, numel(order)));

    T = table(order(:), patchMaxAbsCurv(order), ...
        patchMeanAbsCurv(order), patchMinJacRatio(order), ...
        patchMaxDisp(order), ...
        'VariableNames', {'patch_index', 'max_abs_mean_curv', ...
        'mean_abs_mean_curv', 'min_jac_ratio_to_ref', ...
        'max_displacement_to_ref'});
end

function print_patch_summary(S, S_ref, patchIndex, label)
    inds = S.ixyzs(patchIndex):(S.ixyzs(patchIndex+1)-1);
    [srcvals, ~, ~, ~, ~, ~] = extract_arrays(S);
    [srcRef, ~, ~, ~, ~, ~] = extract_arrays(S_ref);
    jac = surface_jacobian_from_srcvals(srcvals);
    jacRef = surface_jacobian_from_srcvals(srcRef);
    jacRatio = jac ./ max(jacRef, realmin);
    dispToRef = sqrt(sum((S.r - S_ref.r).^2, 1)).';
    curv = S.mean_curv(:);
    [~, localWorst] = max(abs(curv(inds)));

    fprintf('\n  %s, patch %d\n', label, patchIndex);
    fprintf('    worst local/global node = %d/%d\n', ...
        localWorst, inds(localWorst));
    fprintf('    H min/mean/max = %.6e %.6e %.6e\n', ...
        min(curv(inds)), mean(curv(inds)), max(curv(inds)));
    fprintf('    |H| max = %.6e\n', max(abs(curv(inds))));
    fprintf('    jac/ref min/mean/max = %.6e %.6e %.6e\n', ...
        min(jacRatio(inds)), mean(jacRatio(inds)), ...
        max(jacRatio(inds)));
    fprintf('    disp to CAD min/mean/max = %.6e %.6e %.6e\n', ...
        min(dispToRef(inds)), mean(dispToRef(inds)), ...
        max(dispToRef(inds)));
end

function jac = surface_jacobian_from_srcvals(srcvals)
    cr = cross(srcvals(4:6, :).', srcvals(7:9, :).').';
    jac = sqrt(sum(cr.^2, 1)).';
end

function vals = patch_max_abs_curvature(S)
    vals = accumarray(S.patch_id(:), abs(S.mean_curv(:)), ...
        [S.npatches 1], @max, nan);
end

function plot_global_scalar(S, vals, figTitle)
    figure('Name', figTitle);
    plot(S, vals, 'EdgeColor', 'none');
    axis equal;
    view(3);
    grid on;
    camlight headlight;
    lighting gouraud;
    colorbar;
    title(figTitle, 'Interpreter', 'none');
    xlabel('x');
    ylabel('y');
    zlabel('z');
end

function plot_patch_overlay(S_cad, S_primary, S_cad_smooth, patchIndex, nplot)
    [tri, xyzCad, ~] = dense_patch_geometry(S_cad, patchIndex, nplot);
    [~, xyzPrimary, ~] = dense_patch_geometry(S_primary, patchIndex, nplot);
    [~, xyzCadSmooth, ~] = dense_patch_geometry(S_cad_smooth, patchIndex, nplot);

    figure('Name', sprintf('patch %d: CAD / primary / CAD-skeleton overlay', ...
        patchIndex));
    hold on;
    hCad = trisurf(tri, xyzCad(1,:), xyzCad(2,:), xyzCad(3,:), ...
        zeros(size(xyzCad, 2), 1), 'FaceColor', [0.75 0.75 0.75], ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none');
    hPrimary = trisurf(tri, xyzPrimary(1,:), xyzPrimary(2,:), ...
        xyzPrimary(3,:), zeros(size(xyzPrimary, 2), 1), ...
        'FaceColor', [0.1 0.35 0.95], 'FaceAlpha', 0.45, ...
        'EdgeColor', 'none');
    hCadSmooth = trisurf(tri, xyzCadSmooth(1,:), xyzCadSmooth(2,:), ...
        xyzCadSmooth(3,:), zeros(size(xyzCadSmooth, 2), 1), ...
        'FaceColor', [0.95 0.15 0.1], 'FaceAlpha', 0.55, ...
        'EdgeColor', 'none');
    axis equal;
    view(3);
    grid on;
    camlight headlight;
    lighting gouraud;
    title(sprintf('patch %d overlay', patchIndex), 'Interpreter', 'none');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    legend([hCad hPrimary hCadSmooth], ...
        {'CAD input', 'primary smoother', 'CAD-point smoother'}, ...
        'Location', 'best');
end

function plot_patch_diagnostics(S, S_ref, patchIndex, nplot, label)
    [srcvals, ~, ~, ~, ~, ~] = extract_arrays(S);
    [srcRef, ~, ~, ~, ~, ~] = extract_arrays(S_ref);
    jac = surface_jacobian_from_srcvals(srcvals);
    jacRef = surface_jacobian_from_srcvals(srcRef);
    jacRatio = jac ./ max(jacRef, realmin);
    dispToRef = sqrt(sum((S.r - S_ref.r).^2, 1)).';
    logCurv = log10(abs(S.mean_curv(:)) + 1e-16);

    [tri, xyz, uv] = dense_patch_geometry(S, patchIndex, nplot);
    logCurvDense = dense_patch_scalar(S, patchIndex, logCurv, uv);
    jacRatioDense = dense_patch_scalar(S, patchIndex, jacRatio, uv);
    dispDense = dense_patch_scalar(S, patchIndex, dispToRef, uv);

    plot_dense_patch_scalar(tri, xyz, logCurvDense, ...
        sprintf('%s patch %d: log10 |H|', label, patchIndex));
    plot_dense_patch_scalar(tri, xyz, jacRatioDense, ...
        sprintf('%s patch %d: jac/ref', label, patchIndex));
    plot_dense_patch_scalar(tri, xyz, dispDense, ...
        sprintf('%s patch %d: displacement to CAD', label, patchIndex));
end

function [tri, xyz, uv] = dense_patch_geometry(S, patchIndex, nplot)
    order = S.norders(patchIndex);
    inds = S.ixyzs(patchIndex):(S.ixyzs(patchIndex+1)-1);
    uvNodes = surfsmooth3d.internal.koorn.rv_nodes(order);
    valsToCoefs = surfsmooth3d.internal.koorn.vals2coefs(order, uvNodes);
    uv = reference_triangle_grid(nplot);
    tri = delaunay(uv(1,:).', uv(2,:).');
    pols = surfsmooth3d.internal.koorn.pols(order, uv);
    rcoefs = S.r(:, inds) * valsToCoefs.';
    xyz = rcoefs * pols;
end

function valsDense = dense_patch_scalar(S, patchIndex, scalarVals, uv)
    order = S.norders(patchIndex);
    inds = S.ixyzs(patchIndex):(S.ixyzs(patchIndex+1)-1);
    uvNodes = surfsmooth3d.internal.koorn.rv_nodes(order);
    valsToCoefs = surfsmooth3d.internal.koorn.vals2coefs(order, uvNodes);
    pols = surfsmooth3d.internal.koorn.pols(order, uv);
    scoefs = scalarVals(inds).' * valsToCoefs.';
    valsDense = scoefs * pols;
    valsDense = valsDense(:);
end

function uv = reference_triangle_grid(nplot)
    u = [];
    v = [];
    for i = 0:nplot
        for j = 0:(nplot - i)
            u(end+1) = i / nplot; %#ok<AGROW>
            v(end+1) = j / nplot; %#ok<AGROW>
        end
    end
    uv = [u; v];
end

function plot_dense_patch_scalar(tri, xyz, vals, figTitle)
    figure('Name', figTitle);
    trisurf(tri, xyz(1,:), xyz(2,:), xyz(3,:), vals(:), ...
        'EdgeColor', 'none');
    axis equal;
    view(3);
    grid on;
    camlight headlight;
    lighting gouraud;
    colorbar;
    title(figTitle, 'Interpreter', 'none');
    xlabel('x');
    ylabel('y');
    zlabel('z');
end

function results = regression_sweep_examples(repoRoot, doAssert)
%REGRESSION_SWEEP_EXAMPLES Run all bundled STEP examples through MATLAB path.
%
% results = regression_sweep_examples(repoRoot, doAssert)
%
% This is intentionally a numerical smoke/regression test, not a proof of
% full Python parity. It catches the large failures caused by wrong RV
% barycentric conventions, wrong Gmsh rigid-profile settings, and missing
% CAD-edge blending.

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin < 2 || isempty(doAssert)
    doAssert = true;
end

addpath(fullfile(repoRoot, 'matlab'));

stepDir = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files');
files = dir(fullfile(stepDir, '*.step'));
[~, order] = sort({files.name});
files = files(order);

opts = struct();
opts.order = 8;
opts.mesh_fraction = 0.10;
opts.refinement_level = 0;
opts.edge_gll_order = 16;
opts.mesh_profile = 'rigid';

names = strings(numel(files), 1);
npatches = zeros(numel(files), 1);
relAreaError = zeros(numel(files), 1);
maxProjectionDistance = zeros(numel(files), 1);
normalFluxNorm = zeros(numel(files), 1);
projectionFailures = zeros(numel(files), 1);

for k = 1:numel(files)
    names(k) = string(files(k).name);
    stepFile = fullfile(files(k).folder, files(k).name);
    mesh = surfsmooth3d.stepmesher.mesh_step(stepFile, opts);
    stats = surfsmooth3d.stepmesher.validate_area_flux(mesh);

    npatches(k) = mesh.npatches;
    relAreaError(k) = stats.relative_area_error;
    maxProjectionDistance(k) = mesh.meta.max_projection_distance;
    normalFluxNorm(k) = norm(stats.normal_flux);
    projectionFailures(k) = mesh.meta.projection_failures;

    fprintf('%-28s patches=%5d rel=%10.3e fluxnorm=%10.3e maxproj=%10.3e failures=%d\n', ...
        files(k).name, npatches(k), relAreaError(k), normalFluxNorm(k), ...
        maxProjectionDistance(k), projectionFailures(k));
end

results = table(names, npatches, relAreaError, maxProjectionDistance, ...
    normalFluxNorm, projectionFailures);

if doAssert
    assert(all(projectionFailures == 0), 'Projection fallback/failures occurred.');

    torusMask = names == "torus_ellipse.step";
    nonTorus = ~torusMask;

    assert(all(relAreaError(nonTorus) < 2e-3), ...
        'A non-torus bundled example exceeded the regression area tolerance.');

    % The Python reference also reports a large error for this case at
    % order=8, mesh_fraction=0.10, refinement_level=0. Keep it visible, but
    % do not let it silently return to the pre-fix order-one failures.
    assert(all(relAreaError(torusMask) < 0.6), ...
        'torus_ellipse exceeded its known coarse-mesh regression tolerance.');
end
end

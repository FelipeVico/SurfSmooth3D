%% Geometry Package Viewer
% Read-only viewer for step_mesher_cpp_matlab geometry package folders.
%
% A package folder contains:
%   surface.go3
%   scaffold.gidmsh
%   edges_gll.txt
%   metadata.txt

clearvars
close all
clc

%% Locate repository and set paths
driverFile = mfilename('fullpath');
examplesDir = fileparts(driverFile);
repoRoot = fileparts(fileparts(examplesDir));
addpath(fullfile(repoRoot, 'matlab'));

setup_surfsmooth3d();

requiredFunctions = {'surfsmooth3d.surfer.load_from_file', 'surfsmooth3d.internal.koorn.rv_nodes', ...
    'surfsmooth3d.internal.koorn.vals2coefs', 'surfsmooth3d.internal.koorn.pols'};
missingFunctions = requiredFunctions(cellfun(@(f) isempty(which(f)), ...
    requiredFunctions));
if ~isempty(missingFunctions)
    error('stepmesher:missingFmm3dbiePath', ...
        ['This viewer needs SurfSmooth3D/matlab on the MATLAB path. ', ...
         'Missing function(s): %s.'], strjoin(missingFunctions, ', '));
end

%% User parameters
% Geometry family
% packageFamily = 'convergence_test_geometry_manas';
% packageFamily = 'convergence_test_geometry_1';
% packageFamily = 'convergence_test_geometry_2';
% packageFamily = 'convergence_UV_curved_bezier';
packageFamily = 'convergence_test_curved_uv';

% Package stem
packageStem = 'p4_mf0p10000_rl0';
% packageStem = 'p4_mf0p10000_rl1';
% packageStem = 'p4_mf0p10000_rl2';
% packageStem = 'p4_mf0p10000_rl3';

% packageStem = 'p6_mf0p10000_rl0';
% packageStem = 'p6_mf0p10000_rl1';
% packageStem = 'p6_mf0p10000_rl2';
% packageStem = 'p6_mf0p10000_rl3';

% packageStem = 'p8_mf0p10000_rl0';
% packageStem = 'p8_mf0p10000_rl1';
% packageStem = 'p8_mf0p10000_rl2';
% packageStem = 'p8_mf0p10000_rl3';

packageDir = fullfile(repoRoot, 'outputs', 'step', packageFamily, packageStem);


viewerOpts = struct();
viewerOpts.surfaceRefine = 8;
viewerOpts.edgeRefine = 80;
viewerOpts.showSurface = true;
viewerOpts.showEdges = true;
viewerOpts.showEdgeLabels = true;
viewerOpts.showEdgeNodes = false;

%% Load and view
if ~isfolder(packageDir)
    error('stepmesher:fileNotFound', ...
        'Geometry package folder not found: %s', packageDir);
end

pkg = surfsmooth3d.stepmesher.load_geom_package(packageDir);
fprintf('Geometry package: %s\n', packageDir);
fprintf('  surface: %s\n', pkg.surface_file);
fprintf('  scaffold: %s\n', pkg.scaffold_file);
fprintf('  edges: %s\n', pkg.edges_file);
if isfield(pkg, 'edges_format_version') && pkg.edges_format_version < 2
    fprintf('  V1 edge-panel records: %d\n', numel(pkg.edges));
    warning('stepmesher:oldEdgesFile', ...
        ['This package has a V1 edges_gll.txt file, where each record ', ...
         'is one GLL panel rather than one CAD edge. Regenerate the ', ...
         'package to get grouped CAD edges.']);
else
    fprintf('  CAD edges: %d\n', numel(pkg.edges));
end

surfsmooth3d.stepmesher.plot_geom_package(pkg, viewerOpts);

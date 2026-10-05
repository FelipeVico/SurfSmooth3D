%% STEP File Viewer
% Read-only MATLAB viewer for STEP files.
%
% This driver uses MATLAB's Partial Differential Equation Toolbox to import
% a STEP file, display the CAD geometry in 3D, and show face/edge labels.

clearvars
close all
clc

%% Runtime checks
requiredFunctions = {'fegeometry', 'pdegplot'};
missingFunctions = requiredFunctions(cellfun(@(f) exist(f, 'file') ~= 2, ...
    requiredFunctions));

if ~isempty(missingFunctions)
    error('stepmesher:missingToolbox', ...
        ['This viewer requires MATLAB Partial Differential Equation ', ...
         'Toolbox. Missing function(s): %s.'], ...
        strjoin(missingFunctions, ', '));
end

%% Locate repository
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));

%% User parameters
% Choose one STEP file by uncommenting the desired line.
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_geometry_1.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_geometry_2.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_geometry_manas.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'UV_curved_bezier.step');
stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_curved_uv.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'cylinder_bezier_v2.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'intersection_two_balls.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'piecewise_bezier_cap.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'spheres_intersect.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'torus_ellipse.step');
% stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'Multi_res.step');




showFaceLabels = true;
showEdgeLabels = true;
faceAlpha = 0.25;
centerGeometry = true;
axisPadding = 0.08;

% Visual edge mode only. This is not yet an edge-selection tool.
% Options:
%   'labels' - show edge labels when showEdgeLabels is true.
%   'plain'  - plot edges without edge labels.
edgeMode = 'labels';

%% Import and plot
if ~isfile(stepFile)
    error('stepmesher:fileNotFound', 'STEP file not found: %s', stepFile);
end

switch lower(edgeMode)
    case 'labels'
        edgeLabels = on_off(showEdgeLabels);
    case 'plain'
        edgeLabels = "off";
    otherwise
        error('stepmesher:badEdgeMode', ...
            'Unsupported edgeMode "%s". Use "labels" or "plain".', ...
            edgeMode);
end

faceLabels = on_off(showFaceLabels);

fprintf('STEP file: %s\n', stepFile);
fprintf('Importing geometry with fegeometry...\n');
gm = fegeometry(stepFile);
assert(~isempty(gm), 'stepmesher:emptyGeometry', ...
    'fegeometry returned an empty geometry.');

figure('Name', ['STEP viewer: ', get_filename(stepFile)]);
pdegplot(gm, ...
    FaceLabels=faceLabels, ...
    EdgeLabels=edgeLabels, ...
    FaceAlpha=faceAlpha);

view(3)
axis equal
if centerGeometry
    center_geometry_axes(gca, axisPadding);
end
grid on
rotate3d on
title(['STEP viewer: ', get_filename(stepFile)], 'Interpreter', 'none');
xlabel('x');
ylabel('y');
zlabel('z');

fprintf('Viewer ready. Use the mouse to rotate the 3D geometry.\n');

%% Local helpers
function value = on_off(flag)
if flag
    value = "on";
else
    value = "off";
end
end

function name = get_filename(pathname)
[~, name, ext] = fileparts(pathname);
name = [name, ext];
end

function center_geometry_axes(ax, padding)
%CENTER_GEOMETRY_AXES Put plotted 3D geometry in a centered cubic box.
%
% pdegplot can leave the object visually off-center depending on labels,
% aspect ratio, and MATLAB's automatic axis limits. This routine looks at
% the actual plotted XYZ data and sets symmetric limits around its centroid.

if nargin < 2 || isempty(padding)
    padding = 0.08;
end

objects = findall(ax, '-property', 'XData');
xyz = zeros(3, 0);

for i = 1:numel(objects)
    if ~isprop(objects(i), 'YData') || ~isprop(objects(i), 'ZData')
        continue
    end

    x = double(objects(i).XData(:).');
    y = double(objects(i).YData(:).');
    z = double(objects(i).ZData(:).');
    n = min([numel(x), numel(y), numel(z)]);
    if n == 0
        continue
    end

    pts = [x(1:n); y(1:n); z(1:n)];
    pts = pts(:, all(isfinite(pts), 1));
    xyz = [xyz, pts]; %#ok<AGROW>
end

if isempty(xyz)
    axis(ax, 'equal');
    return
end

bmin = min(xyz, [], 2);
bmax = max(xyz, [], 2);
ctr = 0.5 * (bmin + bmax);
span = max(bmax - bmin);
if span <= 0
    span = 1;
end
halfWidth = 0.5 * span * (1 + 2 * padding);

xlim(ax, ctr(1) + halfWidth * [-1, 1]);
ylim(ax, ctr(2) + halfWidth * [-1, 1]);
zlim(ax, ctr(3) + halfWidth * [-1, 1]);
daspect(ax, [1, 1, 1]);
pbaspect(ax, [1, 1, 1]);
axis(ax, 'vis3d');
camtarget(ax, ctr.');
end

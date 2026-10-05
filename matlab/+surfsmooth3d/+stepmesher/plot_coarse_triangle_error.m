function h = plot_coarse_triangle_error(mesh, values, varargin)
%PLOT_COARSE_TRIANGLE_ERROR Color scaffold triangles by one scalar per cell.
%
% h = surfsmooth3d.stepmesher.plot_coarse_triangle_error(mesh, values)
%
% values must contain one entry per coarse scaffold triangle. By default the
% color scale is log10(values), with the colorbar labels shown in the
% original value scale.

p = inputParser;
addParameter(p, 'Title', 'Coarse-triangle error', @(s) ischar(s) || isstring(s));
addParameter(p, 'ColorbarLabel', 'error', @(s) ischar(s) || isstring(s));
addParameter(p, 'Scale', 'log10', @(s) ischar(s) || isstring(s));
addParameter(p, 'Parent', [], @(x) isempty(x) || isa(x, 'matlab.graphics.axis.Axes'));
parse(p, varargin{:});

if ~isfield(mesh, 'scaffold_nodes') || ~isfield(mesh, 'scaffold_triangles')
    error('stepmesher:plotCoarseTriangleError', ...
        'mesh must contain scaffold_nodes and scaffold_triangles.');
end

nodes = double(mesh.scaffold_nodes);
tris = double(mesh.scaffold_triangles);
values = double(values(:));

if size(nodes, 1) ~= 3
    error('stepmesher:plotCoarseTriangleError', ...
        'mesh.scaffold_nodes must have size 3 x nnodes.');
end
if size(tris, 1) ~= 3
    error('stepmesher:plotCoarseTriangleError', ...
        'mesh.scaffold_triangles must have size 3 x ntriangles.');
end
if numel(values) ~= size(tris, 2)
    error('stepmesher:plotCoarseTriangleError', ...
        'values must have one entry per scaffold triangle.');
end

ax = p.Results.Parent;
if isempty(ax)
    ax = gca;
end

[plotValues, colorLimits, tickValues, tickLabels] = scale_face_values( ...
    values, string(p.Results.Scale));

h = patch(ax, ...
    'Vertices', nodes.', ...
    'Faces', tris.', ...
    'FaceVertexCData', plotValues(:), ...
    'FaceColor', 'flat', ...
    'EdgeColor', [0.10, 0.10, 0.10], ...
    'LineWidth', 0.35);

axis(ax, 'equal');
axis(ax, 'vis3d');
grid(ax, 'on');
view(ax, 3);
colormap(ax, turbo(256));
clim(ax, colorLimits);

cb = colorbar(ax);
cb.Label.String = char(p.Results.ColorbarLabel);
if ~isempty(tickValues)
    cb.Ticks = tickValues;
    cb.TickLabels = tickLabels;
end

title(ax, char(p.Results.Title));
xlabel(ax, 'x');
ylabel(ax, 'y');
zlabel(ax, 'z');
rotate3d(ancestor(ax, 'figure'), 'on');
end

function [plotValues, climVals, tickValues, tickLabels] = scale_face_values(values, scaleName)
switch lower(scaleName)
    case "log10"
        finitePositive = values(isfinite(values) & values > 0);
        if isempty(finitePositive)
            floorExp = -16;
            safeValues = 10.^floorExp * ones(size(values));
        else
            floorExp = floor(log10(min(finitePositive))) - 1;
            floorValue = 10.^floorExp;
            safeValues = values;
            safeValues(~isfinite(safeValues) | safeValues <= 0) = floorValue;
        end

        plotValues = log10(safeValues);
        cmin = floor(min(plotValues));
        cmax = ceil(max(plotValues));
        if ~isfinite(cmin)
            cmin = floorExp;
        end
        if ~isfinite(cmax) || cmax <= cmin
            cmax = cmin + 1;
        end
        climVals = [cmin, cmax];

        tickValues = cmin:cmax;
        maxTicks = 8;
        if numel(tickValues) > maxTicks
            tickValues = linspace(cmin, cmax, maxTicks);
        end
        tickLabels = arrayfun(@format_log_tick, tickValues, ...
            'UniformOutput', false);

    otherwise
        finiteValues = values(isfinite(values));
        if isempty(finiteValues)
            plotValues = zeros(size(values));
            climVals = [0, 1];
        else
            fillValue = min(finiteValues);
            plotValues = values;
            plotValues(~isfinite(plotValues)) = fillValue;
            climVals = [min(plotValues), max(plotValues)];
            if climVals(2) <= climVals(1)
                climVals(2) = climVals(1) + 1;
            end
        end
        tickValues = [];
        tickLabels = {};
end
end

function label = format_log_tick(exponent)
if abs(exponent - round(exponent)) < 1e-10
    label = sprintf('1e%d', round(exponent));
else
    label = sprintf('1e%.1f', exponent);
end
end

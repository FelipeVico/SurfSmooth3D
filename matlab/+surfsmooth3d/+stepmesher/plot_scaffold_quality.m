function h = plot_scaffold_quality(mesh, metric, varargin)
%PLOT_SCAFFOLD_QUALITY Color scaffold triangles by a quality diagnostic.
%
% metric can be:
%   'min_angle'  - minimum triangle angle in degrees
%   'edge_ratio' - max edge length / min edge length
%   'area'       - flat scaffold triangle area

if nargin < 2 || isempty(metric)
    metric = 'min_angle';
end
metric = char(metric);

switch metric
    case {'min_angle', 'min_angle_deg'}
        values = mesh.scaffold_min_angle_deg;
        titleText = 'Scaffold minimum angle (degrees)';
    case {'edge_ratio', 'ratio'}
        values = mesh.scaffold_edge_ratio;
        titleText = 'Scaffold max/min edge ratio';
    case 'area'
        values = mesh.scaffold_triangle_area;
        titleText = 'Scaffold triangle area';
    otherwise
        error('stepmesher:plotScaffoldQuality:badMetric', ...
            'Unknown scaffold quality metric: %s', metric);
end

h = surfsmooth3d.stepmesher.plot_coarse_triangle_error(mesh, values, varargin{:});
title(titleText, 'Interpreter', 'none');
end

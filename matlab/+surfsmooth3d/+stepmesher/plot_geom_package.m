function [handles, plotData] = plot_geom_package(pkg, opts)
%PLOT_GEOM_PACKAGE View a STEP mesher geometry package.
%
%   surfsmooth3d.stepmesher.plot_geom_package(pkg) plots the high-order .go3 surface
%   and overlays the CAD edges from edges_gll.txt. Each CAD edge may be a
%   piecewise curve with several GLL panels. The input can be either a
%   package struct returned by surfsmooth3d.stepmesher.load_geom_package or a package
%   directory.

if nargin < 2 || isempty(opts)
    opts = struct();
end

if ischar(pkg) || isstring(pkg)
    pkg = surfsmooth3d.stepmesher.load_geom_package(char(pkg));
end

opts = default_options(opts);
validate_package(pkg);

if ~isfield(pkg, 'surface') || isempty(pkg.surface)
    if ~isempty(which('surfsmooth3d.surfer.load_from_file')) && isfile(pkg.surface_file)
        pkg.surface = surfsmooth3d.surfer.load_from_file(pkg.surface_file);
    else
        error('stepmesher:missingSurfaceLoader', ...
            ['Could not load surface.go3. Add SurfSmooth3D/matlab to the ', ...
             'MATLAB path so surfsmooth3d.surfer.load_from_file is available.']);
    end
end

plotData = struct();
[surfaceVertices, surfaceFaces, patchIds] = build_surface_mesh( ...
    pkg.surface, opts.surfaceRefine);
[edgeCurves, edgeNodeXyz, edgeIds] = build_edge_curves( ...
    pkg.edges, opts.edgeRefine);

plotData.surface_vertices = surfaceVertices;
plotData.surface_faces = surfaceFaces;
plotData.surface_patch_ids = patchIds;
plotData.edge_curves = edgeCurves;
plotData.edge_node_xyz = edgeNodeXyz;
plotData.edge_ids = edgeIds;

fig = figure('Name', figure_name(pkg));
ax = axes('Parent', fig);
hold(ax, 'on');

surfaceVisible = on_off(opts.showSurface);
edgeVisible = on_off(opts.showEdges);
labelVisible = on_off(opts.showEdges && opts.showEdgeLabels);
nodeVisible = on_off(opts.showEdges && opts.showEdgeNodes);

handles = struct();
handles.figure = fig;
handles.axes = ax;
handles.selected_edge_ids = [];

handles.surface = patch(ax, ...
    'Vertices', surfaceVertices, ...
    'Faces', surfaceFaces, ...
    'FaceColor', opts.surfaceColor, ...
    'FaceAlpha', opts.surfaceAlpha, ...
    'EdgeColor', opts.surfaceMeshColor, ...
    'LineWidth', opts.surfaceMeshLineWidth, ...
    'Visible', surfaceVisible, ...
    'HitTest', 'off');

handles.edges = gobjects(1, numel(edgeCurves));
handles.edge_labels = gobjects(1, numel(edgeCurves));
for k = 1:numel(edgeCurves)
    xyz = edgeCurves{k};
    edgeId = edgeIds(k);
    handles.edges(k) = line(ax, xyz(1, :), xyz(2, :), xyz(3, :), ...
        'Color', opts.edgeColor, ...
        'LineWidth', opts.edgeLineWidth, ...
        'Visible', edgeVisible, ...
        'HitTest', 'on', ...
        'PickableParts', 'all', ...
        'ButtonDownFcn', @(src, event) click_edge(src));
    setappdata(handles.edges(k), 'edge_id', edgeId);

    labelXyz = edge_label_point(xyz);
    handles.edge_labels(k) = text(ax, labelXyz(1), labelXyz(2), ...
        labelXyz(3), sprintf('E%d', edgeId), ...
        'FontSize', opts.labelFontSize, ...
        'FontWeight', 'bold', ...
        'Color', opts.labelColor, ...
        'BackgroundColor', opts.labelBackgroundColor, ...
        'Margin', 1, ...
        'Visible', labelVisible, ...
        'HitTest', 'off');
end

if isempty(edgeNodeXyz)
    handles.edge_nodes = gobjects(0);
else
    handles.edge_nodes = scatter3(ax, edgeNodeXyz(1, :), ...
        edgeNodeXyz(2, :), edgeNodeXyz(3, :), ...
        opts.edgeNodeMarkerSize, opts.edgeNodeColor, ...
        'filled', ...
        'Visible', nodeVisible, ...
        'HitTest', 'off');
end

axis(ax, 'equal');
axis(ax, 'vis3d');
grid(ax, 'on');
view(ax, 3);
rotate3d(fig, 'on');
xlabel(ax, 'x');
ylabel(ax, 'y');
zlabel(ax, 'z');
title(ax, figure_name(pkg), 'Interpreter', 'none');
center_axes(ax, surfaceVertices, edgeNodeXyz, opts.axisPadding);

handles.controls = add_controls(fig);
if isfield(pkg, 'edges_format_version') && pkg.edges_format_version < 2
    update_status(sprintf(['Loaded %d V1 edge panels. Regenerate the ', ...
        'package for true CAD-edge grouping.'], numel(edgeCurves)));
else
    update_status(sprintf('Loaded %d CAD edges. Click an edge or enter IDs.', ...
        numel(edgeCurves)));
end

    function controls = add_controls(parentFig)
        panel = uipanel('Parent', parentFig, ...
            'Title', 'Geometry package viewer', ...
            'Units', 'normalized', ...
            'Position', [0.015, 0.015, 0.30, 0.25]);

        controls.surface = uicontrol(panel, ...
            'Style', 'checkbox', ...
            'String', 'Surface', ...
            'Units', 'normalized', ...
            'Position', [0.04, 0.78, 0.28, 0.16], ...
            'Value', logical(opts.showSurface), ...
            'Callback', @(src, event) toggle_surface(src));

        controls.edges = uicontrol(panel, ...
            'Style', 'checkbox', ...
            'String', 'Edges', ...
            'Units', 'normalized', ...
            'Position', [0.34, 0.78, 0.25, 0.16], ...
            'Value', logical(opts.showEdges), ...
            'Callback', @(src, event) toggle_edges(src));

        controls.labels = uicontrol(panel, ...
            'Style', 'checkbox', ...
            'String', 'Labels', ...
            'Units', 'normalized', ...
            'Position', [0.62, 0.78, 0.30, 0.16], ...
            'Value', logical(opts.showEdgeLabels), ...
            'Callback', @(src, event) toggle_labels(src));

        controls.nodes = uicontrol(panel, ...
            'Style', 'checkbox', ...
            'String', 'GLL nodes', ...
            'Units', 'normalized', ...
            'Position', [0.04, 0.58, 0.40, 0.16], ...
            'Value', logical(opts.showEdgeNodes), ...
            'Callback', @(src, event) toggle_nodes(src));

        uicontrol(panel, ...
            'Style', 'text', ...
            'String', 'Edge IDs', ...
            'HorizontalAlignment', 'left', ...
            'Units', 'normalized', ...
            'Position', [0.04, 0.38, 0.24, 0.14]);

        controls.edgeEdit = uicontrol(panel, ...
            'Style', 'edit', ...
            'String', '', ...
            'HorizontalAlignment', 'left', ...
            'Units', 'normalized', ...
            'Position', [0.28, 0.40, 0.66, 0.14]);

        controls.highlight = uicontrol(panel, ...
            'Style', 'pushbutton', ...
            'String', 'Highlight', ...
            'Units', 'normalized', ...
            'Position', [0.04, 0.20, 0.42, 0.14], ...
            'Callback', @(src, event) highlight_from_edit());

        controls.clear = uicontrol(panel, ...
            'Style', 'pushbutton', ...
            'String', 'Clear', ...
            'Units', 'normalized', ...
            'Position', [0.52, 0.20, 0.42, 0.14], ...
            'Callback', @(src, event) set_selected_edges([]));

        controls.status = uicontrol(panel, ...
            'Style', 'text', ...
            'String', '', ...
            'HorizontalAlignment', 'left', ...
            'Units', 'normalized', ...
            'Position', [0.04, 0.02, 0.90, 0.14]);
    end

    function toggle_surface(src)
        set(handles.surface, 'Visible', on_off(get(src, 'Value')));
    end

    function toggle_edges(src)
        value = logical(get(src, 'Value'));
        set(handles.edges, 'Visible', on_off(value));
        if ~isempty(handles.edge_nodes) && all(isgraphics(handles.edge_nodes))
            set(handles.edge_nodes, 'Visible', on_off(value && ...
                logical(get(handles.controls.nodes, 'Value'))));
        end
        set(handles.edge_labels, 'Visible', on_off(value && ...
            logical(get(handles.controls.labels, 'Value'))));
    end

    function toggle_labels(src)
        value = logical(get(src, 'Value')) && ...
            logical(get(handles.controls.edges, 'Value'));
        set(handles.edge_labels, 'Visible', on_off(value));
    end

    function toggle_nodes(src)
        if ~isempty(handles.edge_nodes) && all(isgraphics(handles.edge_nodes))
            value = logical(get(src, 'Value')) && ...
                logical(get(handles.controls.edges, 'Value'));
            set(handles.edge_nodes, 'Visible', on_off(value));
        end
    end

    function highlight_from_edit()
        ids = parse_edge_ids(get(handles.controls.edgeEdit, 'String'), ...
            edgeIds);
        set_selected_edges(ids);
    end

    function click_edge(src)
        edgeId = getappdata(src, 'edge_id');
        set(handles.controls.edgeEdit, 'String', sprintf('%d', edgeId));
        set_selected_edges(edgeId);
    end

    function set_selected_edges(ids)
        handles.selected_edge_ids = unique(ids(:).');
        reset_edge_styles();

        if isempty(handles.selected_edge_ids)
            update_status('No selected edges.');
            return
        end

        found = false(size(handles.selected_edge_ids));
        for ii = 1:numel(handles.selected_edge_ids)
            idx = find(edgeIds == handles.selected_edge_ids(ii), 1);
            if isempty(idx)
                continue
            end
            found(ii) = true;
            set(handles.edges(idx), ...
                'Color', opts.highlightEdgeColor, ...
                'LineWidth', opts.highlightEdgeLineWidth);
            set(handles.edge_labels(idx), ...
                'Color', opts.highlightEdgeColor, ...
                'FontWeight', 'bold');
        end

        if all(found)
            update_status(sprintf('Selected edge(s): %s', ...
                join_ids(handles.selected_edge_ids)));
        else
            missing = handles.selected_edge_ids(~found);
            update_status(sprintf('Selected %s. Missing ID(s): %s', ...
                join_ids(handles.selected_edge_ids(found)), ...
                join_ids(missing)));
        end
    end

    function reset_edge_styles()
        for ii = 1:numel(handles.edges)
            set(handles.edges(ii), ...
                'Color', opts.edgeColor, ...
                'LineWidth', opts.edgeLineWidth);
            set(handles.edge_labels(ii), ...
                'Color', opts.labelColor, ...
                'FontWeight', 'bold');
        end
    end

    function update_status(message)
        if isfield(handles, 'controls') && isfield(handles.controls, 'status')
            set(handles.controls.status, 'String', message);
        end
        fprintf('%s\n', message);
    end
end

function opts = default_options(opts)
defaults = struct();
defaults.surfaceRefine = 8;
defaults.edgeRefine = 80;
defaults.showSurface = true;
defaults.showEdges = true;
defaults.showEdgeLabels = true;
defaults.showEdgeNodes = false;
defaults.surfaceColor = [0.72, 0.78, 0.84];
defaults.surfaceAlpha = 0.88;
defaults.surfaceMeshColor = [0.58, 0.62, 0.66];
defaults.surfaceMeshLineWidth = 0.25;
defaults.edgeColor = [0.02, 0.02, 0.02];
defaults.edgeLineWidth = 1.1;
defaults.highlightEdgeColor = [0.92, 0.12, 0.06];
defaults.highlightEdgeLineWidth = 3.0;
defaults.labelColor = [0.00, 0.00, 0.00];
defaults.labelBackgroundColor = [1.00, 1.00, 1.00];
defaults.labelFontSize = 9;
defaults.edgeNodeColor = [0.95, 0.18, 0.10];
defaults.edgeNodeMarkerSize = 16;
defaults.axisPadding = 0.08;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end

function validate_package(pkg)
required = {'package_dir', 'surface_file', 'edges_file', 'edges'};
for i = 1:numel(required)
    if ~isfield(pkg, required{i})
        error('stepmesher:badPackage', ...
            'Package is missing field "%s".', required{i});
    end
end
if ~isfile(pkg.surface_file)
    error('stepmesher:fileNotFound', ...
        'Package surface file not found: %s', pkg.surface_file);
end
if ~isfile(pkg.edges_file)
    error('stepmesher:fileNotFound', ...
        'Package edge file not found: %s', pkg.edges_file);
end
if isempty(which('surfsmooth3d.internal.koorn.rv_nodes')) || isempty(which('surfsmooth3d.internal.koorn.pols'))
    error('stepmesher:missingKoorn', ...
        'Add SurfSmooth3D/matlab to the MATLAB path so koorn is available.');
end
end

function [vertices, faces, patchIds] = build_surface_mesh(S, refine)
[srcvals, ~, norders, ixyzs, iptype] = extract_arrays(S);
if any(iptype ~= 1)
    error('stepmesher:unsupportedPatchType', ...
        'plot_geom_package currently supports triangular RV patches only.');
end

refine = max(1, round(double(refine)));
[uvPlot, triFaces] = reference_triangle_mesh(refine);
npatches = numel(norders);
nplot = size(uvPlot, 2);
nfaces = size(triFaces, 1);

vertices = zeros(npatches * nplot, 3);
faces = zeros(npatches * nfaces, 3);
patchIds = zeros(npatches * nfaces, 1);

cache = struct();
for patch = 1:npatches
    order = double(norders(patch));
    cacheKey = sprintf('o%d', order);
    if ~isfield(cache, cacheKey)
        uvRV = surfsmooth3d.internal.koorn.rv_nodes(order);
        cache.(cacheKey).A = surfsmooth3d.internal.koorn.vals2coefs(order, uvRV);
        cache.(cacheKey).P = surfsmooth3d.internal.koorn.pols(order, uvPlot);
    end

    nodeCols = ixyzs(patch):(ixyzs(patch + 1) - 1);
    R = srcvals(1:3, nodeCols);
    coef = R * cache.(cacheKey).A.';
    rPlot = coef * cache.(cacheKey).P;

    vRange = (patch - 1) * nplot + (1:nplot);
    fRange = (patch - 1) * nfaces + (1:nfaces);
    vertices(vRange, :) = rPlot.';
    faces(fRange, :) = triFaces + vRange(1) - 1;
    patchIds(fRange) = patch;
end

bad = ~all(isfinite(vertices), 2);
if any(bad)
    error('stepmesher:badSurface', ...
        'Dense surface evaluation produced nonfinite vertices.');
end
end

function [uv, faces] = reference_triangle_mesh(refine)
index = nan(refine + 1, refine + 1);
uv = zeros(2, (refine + 1) * (refine + 2) / 2);

k = 0;
for i = 0:refine
    for j = 0:(refine - i)
        k = k + 1;
        index(i + 1, j + 1) = k;
        uv(:, k) = [i; j] / refine;
    end
end

faces = zeros(refine * refine, 3);
nf = 0;
for i = 0:(refine - 1)
    for j = 0:(refine - 1 - i)
        v00 = index(i + 1, j + 1);
        v10 = index(i + 2, j + 1);
        v01 = index(i + 1, j + 2);
        nf = nf + 1;
        faces(nf, :) = [v00, v10, v01];

        if i + j <= refine - 2
            v11 = index(i + 2, j + 2);
            nf = nf + 1;
            faces(nf, :) = [v10, v11, v01];
        end
    end
end
faces = faces(1:nf, :);
end

function [edgeCurves, edgeNodeXyz, edgeIds] = build_edge_curves(edges, refine)
nedges = numel(edges);
edgeCurves = cell(1, nedges);
edgeIds = zeros(1, nedges);
edgeNodeXyz = zeros(3, 0);

for k = 1:nedges
    edge = edges(k);
    edgeIds(k) = edge.edge_id;
    curve = zeros(3, 0);
    for panelIdx = 1:edge.num_panels
        panel = edge.panels(panelIdx);
        nplot = max([panel.num_nodes, round(double(refine)), 2]);
        aPlot = linspace(min(panel.a), max(panel.a), nplot);
        panelCurve = barycentric_interp_1d(panel.a, panel.xyz, aPlot);
        if isempty(curve)
            curve = panelCurve;
        else
            curve = [curve, nan(3, 1), panelCurve]; %#ok<AGROW>
        end
        edgeNodeXyz = [edgeNodeXyz, panel.xyz]; %#ok<AGROW>
    end
    edgeCurves{k} = curve;
end
end

function xyz = edge_label_point(curve)
finiteCols = find(all(isfinite(curve), 1));
if isempty(finiteCols)
    xyz = [0; 0; 0];
else
    xyz = curve(:, finiteCols(max(1, round(numel(finiteCols) / 2))));
end
end

function valuesPlot = barycentric_interp_1d(nodes, values, query)
nodes = double(nodes(:));
query = double(query(:).');
values = double(values);
n = numel(nodes);

if size(values, 2) ~= n
    error('stepmesher:badEdges', ...
        'Number of interpolation nodes does not match edge xyz data.');
end

w = ones(n, 1);
for j = 1:n
    others = [1:(j - 1), (j + 1):n];
    w(j) = 1.0 / prod(nodes(j) - nodes(others));
end

valuesPlot = zeros(size(values, 1), numel(query));
tol = 100 * eps(max(1, max(abs(nodes))));
for iq = 1:numel(query)
    diff = query(iq) - nodes;
    hit = find(abs(diff) <= tol, 1);
    if ~isempty(hit)
        valuesPlot(:, iq) = values(:, hit);
    else
        weights = w ./ diff;
        valuesPlot(:, iq) = values * (weights / sum(weights));
    end
end
end

function center_axes(ax, surfaceVertices, edgeNodeXyz, padding)
xyz = surfaceVertices.';
if ~isempty(edgeNodeXyz)
    xyz = [xyz, edgeNodeXyz];
end
xyz = xyz(:, all(isfinite(xyz), 1));
if isempty(xyz)
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
camtarget(ax, ctr.');
end

function value = on_off(flag)
if flag
    value = 'on';
else
    value = 'off';
end
end

function ids = parse_edge_ids(textValue, allIds)
textValue = strtrim(char(textValue));
if isempty(textValue)
    ids = [];
    return
end
if strcmpi(textValue, 'all')
    ids = allIds;
    return
end

clean = regexprep(textValue, '[\[\],;]', ' ');
ids = sscanf(clean, '%f').';
ids = unique(round(ids(isfinite(ids))));
end

function textValue = join_ids(ids)
if isempty(ids)
    textValue = 'none';
else
    textValue = char(strjoin(compose('%d', ids), ', '));
end
end

function name = figure_name(pkg)
[~, folderName] = fileparts(pkg.package_dir);
name = ['Geometry package: ', folderName];
end

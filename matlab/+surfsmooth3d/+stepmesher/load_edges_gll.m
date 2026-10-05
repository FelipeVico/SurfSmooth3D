function edges = load_edges_gll(edgeFile)
%LOAD_EDGES_GLL Load STEP_MESHER_EDGES_GLL edge text files.
%
% V2 files return one struct per CAD edge, with one or more GLL panels in
% edges(k).panels. V1 panel-level files are accepted as a compatibility
% fallback and normalized to the same shape.

fid = fopen(edgeFile, 'r');
if fid < 0
    error('stepmesher:io', 'Could not open %s for reading.', edgeFile);
end
cleanup = onCleanup(@() fclose(fid));

marker = next_nonempty_line(fid);
switch marker
    case 'STEP_MESHER_EDGES_GLL_V2'
        edges = read_v2(fid, edgeFile);
    case 'STEP_MESHER_EDGES_GLL_V1'
        edges = read_v1_as_v2(fid, edgeFile);
    otherwise
        error('stepmesher:badEdgesFile', ...
            'Unsupported edge file marker in %s.', edgeFile);
end
end

function edges = read_v2(fid, edgeFile)
edgeOrder = read_scalar_line(fid, 'edge_order');
numEdges = read_scalar_line(fid, 'num_cad_edges');
edges = repmat(empty_edge(edgeOrder), 1, numEdges);

for k = 1:numEdges
    line = next_nonempty_line(fid);
    if ~strcmp(line, 'CAD_EDGE')
        error('stepmesher:badEdgesFile', ...
            'Expected CAD_EDGE block %d in %s.', k, edgeFile);
    end

    edges(k).edge_order = edgeOrder;
    edges(k).format_version = 2;
    edges(k).edge_id = read_scalar_line(fid, 'edge_id');
    edges(k).curve_tag = read_scalar_line(fid, 'curve_tag');
    edges(k).adjacent_faces = read_vector_line(fid, 'adjacent_faces', 2);
    edges(k).num_panels = read_scalar_line(fid, 'num_panels');
    edges(k).panels = repmat(empty_panel(edgeOrder), ...
        1, edges(k).num_panels);

    for j = 1:edges(k).num_panels
        edges(k).panels(j) = read_panel(fid, edgeFile, edgeOrder);
    end

    line = next_nonempty_line(fid);
    if ~strcmp(line, 'END_CAD_EDGE')
        error('stepmesher:badEdgesFile', ...
            'Expected END_CAD_EDGE after CAD edge %d.', k);
    end
end

line = next_nonempty_line(fid);
if ~strcmp(line, 'END_STEP_MESHER_EDGES_GLL')
    error('stepmesher:badEdgesFile', ...
        'Expected END_STEP_MESHER_EDGES_GLL in %s.', edgeFile);
end
end

function edges = read_v1_as_v2(fid, edgeFile)
edgeOrder = read_scalar_line(fid, 'edge_order');
numEdges = read_scalar_line(fid, 'num_edges');
edges = repmat(empty_edge(edgeOrder), 1, numEdges);

for k = 1:numEdges
    line = next_nonempty_line(fid);
    if ~strcmp(line, 'EDGE')
        error('stepmesher:badEdgesFile', ...
            'Expected EDGE block %d in %s.', k, edgeFile);
    end

    edges(k).edge_order = edgeOrder;
    edges(k).format_version = 1;
    edges(k).edge_id = read_scalar_line(fid, 'edge_id');
    edges(k).curve_tag = read_scalar_line(fid, 'curve_tag');
    nodeTags = read_vector_line(fid, 'node_tags', 2);
    edges(k).adjacent_faces = read_vector_line(fid, 'adjacent_faces', 2);
    numNodes = read_scalar_line(fid, 'num_nodes');
    columns = next_nonempty_line(fid);
    if ~startsWith(columns, 'columns ')
        error('stepmesher:badEdgesFile', ...
            'Expected columns line for V1 edge %d.', k);
    end

    data = read_panel_data(fid, numNodes, k);
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

    line = next_nonempty_line(fid);
    if ~strcmp(line, 'END_EDGE')
        error('stepmesher:badEdgesFile', ...
            'Expected END_EDGE after V1 edge %d.', k);
    end
end

line = next_nonempty_line(fid);
if ~strcmp(line, 'END_STEP_MESHER_EDGES_GLL')
    error('stepmesher:badEdgesFile', ...
        'Expected END_STEP_MESHER_EDGES_GLL in %s.', edgeFile);
end
end

function panel = read_panel(fid, edgeFile, edgeOrder)
line = next_nonempty_line(fid);
if ~strcmp(line, 'PANEL')
    error('stepmesher:badEdgesFile', ...
        'Expected PANEL block in %s.', edgeFile);
end

panel = empty_panel(edgeOrder);
panel.panel_id = read_scalar_line(fid, 'panel_id');
panel.node_tags = read_vector_line(fid, 'node_tags', 2);
panel.num_nodes = read_scalar_line(fid, 'num_nodes');

columns = next_nonempty_line(fid);
if ~startsWith(columns, 'columns ')
    error('stepmesher:badEdgesFile', ...
        'Expected columns line for panel %d.', panel.panel_id);
end

data = read_panel_data(fid, panel.num_nodes, panel.panel_id);
panel.a = data(:, 1).';
panel.xyz = data(:, 2:4).';
panel.tangent = data(:, 5:7).';
panel.w = data(:, 8).';

line = next_nonempty_line(fid);
if ~strcmp(line, 'END_PANEL')
    error('stepmesher:badEdgesFile', ...
        'Expected END_PANEL after panel %d.', panel.panel_id);
end
end

function data = read_panel_data(fid, numNodes, panelId)
data = fscanf(fid, '%f', [8, numNodes]).';
if size(data, 1) ~= numNodes || size(data, 2) ~= 8
    error('stepmesher:badEdgesFile', ...
        'Could not read %d edge nodes for panel %d.', numNodes, panelId);
end
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

function line = next_nonempty_line(fid)
line = '';
while ischar(line)
    line = fgetl(fid);
    if ~ischar(line)
        error('stepmesher:badEdgesFile', 'Unexpected end of file.');
    end
    line = strtrim(line);
    if ~isempty(line)
        return
    end
end
end

function value = read_scalar_line(fid, key)
values = read_vector_line(fid, key, 1);
value = values(1);
end

function values = read_vector_line(fid, key, n)
line = next_nonempty_line(fid);
parts = strsplit(line);
if numel(parts) ~= n + 1 || ~strcmp(parts{1}, key)
    error('stepmesher:badEdgesFile', 'Expected "%s" line.', key);
end
values = zeros(1, n);
for i = 1:n
    values(i) = str2double(parts{i + 1});
end
end

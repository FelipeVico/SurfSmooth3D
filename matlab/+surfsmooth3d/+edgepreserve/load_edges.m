function edges = load_edges(edgeFile)
    loader = which('surfsmooth3d.stepmesher.load_edges_gll');
    if ~isempty(loader)
        edges = surfsmooth3d.stepmesher.load_edges_gll(edgeFile);
        return
    end

    fid = fopen(edgeFile, 'r');
    if fid < 0
        error('edgepreserve:edgesOpenFailed', ...
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
            error('edgepreserve:badEdgesFile', ...
                'Unsupported edge file marker in %s.', edgeFile);
    end
    clear cleanup;
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
            error('edgepreserve:badEdgesFile', ...
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
        error('edgepreserve:badEdgesFile', ...
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
        error('edgepreserve:badEdgesFile', ...
            'Expected line beginning with "%s".', key);
    end
    vals = str2double(parts(2:end)).';
end


function expect_line(fid, expected, edgeFile)
    line = next_nonempty_line(fid);
    if ~strcmp(line, expected)
        error('edgepreserve:badEdgesFile', ...
            'Expected "%s" in %s, got "%s".', expected, edgeFile, line);
    end
end


function line = next_nonempty_line(fid)
    while true
        line = fgetl(fid);
        if ~ischar(line)
            error('edgepreserve:unexpectedEof', ...
                'Unexpected end of edge file.');
        end
        line = strtrim(line);
        if ~isempty(line)
            return
        end
    end
end

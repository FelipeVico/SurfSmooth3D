function payload = export_edges_gll(mesh, outFile)
%EXPORT_EDGES_GLL Write CAD edges with piecewise GLL panels.
%
% The V2 output stores one top-level record per CAD/Gmsh curve. Each CAD
% edge contains one or more GLL panels corresponding to the 1D mesh chunks
% along that curve.

required = {'edge_order', 'edge_nodes_per_panel', 'ncad_edges', ...
    'npanels_total', 'edge_gll_nodes', 'edge_ids', 'edge_curve_tags', ...
    'edge_adjacent_faces', 'edge_panel_start', 'edge_panel_count', ...
    'panel_node_tags', 'panel_xyz', 'panel_tangents', 'panel_weights'};
for k = 1:numel(required)
    if ~isfield(mesh, required{k})
        error('stepmesher:badMesh', 'Mesh is missing field "%s".', ...
            required{k});
    end
end

edgeOrder = double(mesh.edge_order);
nodesPerPanel = double(mesh.edge_nodes_per_panel);
numCadEdges = double(mesh.ncad_edges);
numPanels = double(mesh.npanels_total);
totalPanelNodes = nodesPerPanel * numPanels;

if edgeOrder + 1 ~= nodesPerPanel
    error('stepmesher:badEdges', ...
        'edge_nodes_per_panel must equal edge_order + 1.');
end
if numCadEdges < 0 || numPanels < 0 || ...
        totalPanelNodes ~= size(mesh.panel_xyz, 2)
    error('stepmesher:badEdges', ...
        'Panel xyz array size is inconsistent with edge metadata.');
end
if size(mesh.panel_xyz, 1) ~= 3 || size(mesh.panel_tangents, 1) ~= 3
    error('stepmesher:badEdges', ...
        'panel_xyz and panel_tangents must have size 3 x total nodes.');
end
if size(mesh.panel_tangents, 2) ~= totalPanelNodes || ...
        numel(mesh.panel_weights) ~= totalPanelNodes
    error('stepmesher:badEdges', ...
        'Panel tangent/weight arrays are inconsistent with panel_xyz.');
end
if numel(mesh.edge_ids) ~= numCadEdges || ...
        numel(mesh.edge_curve_tags) ~= numCadEdges || ...
        size(mesh.edge_adjacent_faces, 1) ~= 2 || ...
        size(mesh.edge_adjacent_faces, 2) ~= numCadEdges || ...
        numel(mesh.edge_panel_start) ~= numCadEdges || ...
        numel(mesh.edge_panel_count) ~= numCadEdges || ...
        size(mesh.panel_node_tags, 1) ~= 2 || ...
        size(mesh.panel_node_tags, 2) ~= numPanels
    error('stepmesher:badEdges', ...
        'CAD edge metadata arrays are inconsistent.');
end

fid = fopen(outFile, 'w');
if fid < 0
    error('stepmesher:io', 'Could not open %s for writing.', outFile);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'STEP_MESHER_EDGES_GLL_V2\n');
fprintf(fid, 'edge_order %d\n', edgeOrder);
fprintf(fid, 'num_cad_edges %d\n\n', numCadEdges);

for edgeIdx = 1:numCadEdges
    panelStart = double(mesh.edge_panel_start(edgeIdx));
    panelCount = double(mesh.edge_panel_count(edgeIdx));
    if panelStart < 1 || panelCount < 1 || ...
            panelStart + panelCount - 1 > numPanels
        error('stepmesher:badEdges', ...
            'Invalid panel range for CAD edge %d.', edgeIdx);
    end

    fprintf(fid, 'CAD_EDGE\n');
    fprintf(fid, 'edge_id %d\n', mesh.edge_ids(edgeIdx));
    fprintf(fid, 'curve_tag %d\n', mesh.edge_curve_tags(edgeIdx));
    fprintf(fid, 'adjacent_faces %d %d\n', ...
        mesh.edge_adjacent_faces(1, edgeIdx), ...
        mesh.edge_adjacent_faces(2, edgeIdx));
    fprintf(fid, 'num_panels %d\n\n', panelCount);

    for localPanel = 1:panelCount
        panelIdx = panelStart + localPanel - 1;
        cols = (panelIdx - 1) * nodesPerPanel + (1:nodesPerPanel);

        fprintf(fid, 'PANEL\n');
        fprintf(fid, 'panel_id %d\n', localPanel);
        fprintf(fid, 'node_tags %d %d\n', ...
            mesh.panel_node_tags(1, panelIdx), ...
            mesh.panel_node_tags(2, panelIdx));
        fprintf(fid, 'num_nodes %d\n', nodesPerPanel);
        fprintf(fid, 'columns a x y z tx ty tz w\n');

        for j = 1:nodesPerPanel
            col = cols(j);
            fprintf(fid, ['%.16e %.16e %.16e %.16e %.16e %.16e %.16e ', ...
                '%.16e\n'], mesh.edge_gll_nodes(j), ...
                mesh.panel_xyz(1, col), mesh.panel_xyz(2, col), ...
                mesh.panel_xyz(3, col), mesh.panel_tangents(1, col), ...
                mesh.panel_tangents(2, col), mesh.panel_tangents(3, col), ...
                mesh.panel_weights(col));
        end
        fprintf(fid, 'END_PANEL\n\n');
    end

    fprintf(fid, 'END_CAD_EDGE\n\n');
end

fprintf(fid, 'END_STEP_MESHER_EDGES_GLL\n');

payload = struct();
payload.edges_file = outFile;
payload.edge_order = edgeOrder;
payload.num_cad_edges = numCadEdges;
payload.num_panels = numPanels;
payload.nodes_per_panel = nodesPerPanel;
payload.total_panel_nodes = totalPanelNodes;
payload.edge_gll_nodes = mesh.edge_gll_nodes;
end

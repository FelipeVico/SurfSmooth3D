function payload = export_gidmsh(mesh, outFile)
%EXPORT_GIDMSH Write the flat scaffold mesh in fmm3dbie GiD .gidmsh format.
%
% The surface smoother reads this as ifiletype = 3. It expects flat
% 3-node triangles; the smoother adds quadratic midpoint nodes internally.

required = {'scaffold_nodes', 'scaffold_triangles'};
for k = 1:numel(required)
    if ~isfield(mesh, required{k})
        error('stepmesher:badMesh', 'Mesh is missing field "%s".', required{k});
    end
end

nodes = double(mesh.scaffold_nodes);
tris = double(mesh.scaffold_triangles);

if size(nodes, 1) ~= 3
    error('stepmesher:badMesh', ...
        'mesh.scaffold_nodes must have size 3 x nnodes.');
end
if size(tris, 1) ~= 3
    error('stepmesher:badMesh', ...
        'mesh.scaffold_triangles must have size 3 x ntriangles.');
end
if any(tris(:) ~= round(tris(:))) || any(tris(:) < 1) || ...
        any(tris(:) > size(nodes, 2))
    error('stepmesher:badMesh', ...
        'mesh.scaffold_triangles must contain valid 1-based node ids.');
end

fid = fopen(outFile, 'w');
if fid < 0
    error('stepmesher:io', 'Could not open %s for writing.', outFile);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'MESH dimension 3 ElemType Triangle Nnode 3\n');
fprintf(fid, 'Coordinates\n');
for i = 1:size(nodes, 2)
    fprintf(fid, '%d %.16e %.16e %.16e\n', i, nodes(1, i), nodes(2, i), ...
        nodes(3, i));
end
fprintf(fid, 'End Coordinates\n');
fprintf(fid, 'Elements\n');
for i = 1:size(tris, 2)
    fprintf(fid, '%d %d %d %d\n', i, tris(1, i), tris(2, i), tris(3, i));
end
fprintf(fid, 'End Elements\n');

payload = struct();
payload.gidmsh_file = outFile;
payload.num_nodes = size(nodes, 2);
payload.num_triangles = size(tris, 2);
end

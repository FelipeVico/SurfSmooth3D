function payload = export_xyzw(mesh, outFile)
%EXPORT_XYZW Write Python-reference compatible x,y,z,w text output.
%
% The first line is "order,num_triangles"; following rows are x,y,z,w with
% comma delimiters.

[srcvals, ~, ~, meta] = surfsmooth3d.stepmesher.to_srcvals(mesh);
order = double(mesh.norder);
npatches = double(mesh.npatches);
nodes_per_patch = double(mesh.nodes_per_patch);
wref = surfsmooth3d.internal.koorn.rv_weights(order);
if numel(wref) ~= nodes_per_patch
    error('stepmesher:badOrder', ...
        'surfsmooth3d.internal.koorn.rv_weights(%d) gives %d weights, but mesh has %d nodes.', ...
        order, numel(wref), nodes_per_patch);
end
wref = wref(:).';

du = srcvals(4:6, :);
dv = srcvals(7:9, :);
jac = vecnorm(cross(du, dv, 1), 2, 1);
weights = zeros(1, size(srcvals, 2));
for patch = 1:npatches
    cols = (patch - 1) * nodes_per_patch + (1:nodes_per_patch);
    weights(cols) = jac(cols) .* wref;
end

xyzw = [srcvals(1:3, :).', weights(:)];

fid = fopen(outFile, 'w');
if fid < 0
    error('stepmesher:io', 'Could not open %s for writing.', outFile);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%d,%d\n', order, npatches);
fprintf(fid, '%.16e,%.16e,%.16e,%.16e\n', xyzw.');

payload = struct();
payload.xyzw = xyzw;
payload.weights = weights(:);
payload.num_triangles = npatches;
payload.nodes_per_tri = nodes_per_patch;
payload.total_nodes = size(srcvals, 2);
payload.meta = meta;
end

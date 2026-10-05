function payload = export_go3(mesh, outFile)
%EXPORT_GO3 Write fmm3dbie .go3 high-order surface output.
%
% The .go3 file stores one common triangular RV order, the number of
% patches, and srcvals = [r; du; dv; n] in fmm3dbie column-major order.

[srcvals, norders, iptype, meta] = surfsmooth3d.stepmesher.to_srcvals(mesh);

if isempty(norders)
    error('stepmesher:badMesh', 'mesh contains no patches.');
end

order = double(norders(1));
npatches = double(numel(norders));
nodes_per_patch = double(mesh.nodes_per_patch);
expected_nodes = (order + 1) * (order + 2) / 2;

if order ~= round(order) || order < 0
    error('stepmesher:badOrder', ...
        '.go3 export requires a nonnegative integer order.');
end
if any(double(norders(:)) ~= order)
    error('stepmesher:badOrder', ...
        '.go3 export requires one common order for all patches.');
end
if any(double(iptype(:)) ~= 1)
    error('stepmesher:badPatchType', ...
        '.go3 export currently supports triangular RV patches only.');
end
if nodes_per_patch ~= expected_nodes
    error('stepmesher:badOrder', ...
        ['.go3 export expected %d nodes per patch for order %d, ', ...
         'but mesh has %d.'], expected_nodes, order, nodes_per_patch);
end

fid = fopen(outFile, 'w');
if fid < 0
    error('stepmesher:io', 'Could not open %s for writing.', outFile);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '%d\n', order);
fprintf(fid, '%d\n', npatches);
fprintf(fid, '%.16e\n', srcvals.');

payload = struct();
payload.go3_file = outFile;
payload.norder = order;
payload.npatches = npatches;
payload.total_nodes = size(srcvals, 2);
payload.srcvals = srcvals;
payload.meta = meta;
end

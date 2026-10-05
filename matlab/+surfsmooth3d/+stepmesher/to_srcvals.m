function [srcvals, norders, iptype, meta] = to_srcvals(mesh)
%TO_SRCVALS Convert projected patch nodes to fmm3dbie srcvals.

required = {'xyz', 'nodes_per_patch', 'npatches', 'norder'};
for k = 1:numel(required)
    if ~isfield(mesh, required{k})
        error('stepmesher:badMesh', 'Mesh is missing field "%s".', required{k});
    end
end

order = double(mesh.norder);
npatches = double(mesh.npatches);
nodes_per_patch = double(mesh.nodes_per_patch);
xyz = double(mesh.xyz);

if size(xyz, 1) ~= 3
    error('stepmesher:badMesh', 'mesh.xyz must have size 3 x total_nodes.');
end
if size(xyz, 2) ~= npatches * nodes_per_patch
    error('stepmesher:badMesh', ...
        'mesh.xyz columns do not match npatches * nodes_per_patch.');
end

uv = surfsmooth3d.internal.koorn.rv_nodes(order);
if size(uv, 2) ~= nodes_per_patch
    error('stepmesher:badOrder', ...
        'surfsmooth3d.internal.koorn.rv_nodes(%d) gives %d nodes, but mesh has %d.', ...
        order, size(uv, 2), nodes_per_patch);
end

A = surfsmooth3d.internal.koorn.vals2coefs(order, uv);
[~, Du, Dv] = surfsmooth3d.internal.koorn.ders(order, uv);

srcvals = zeros(12, size(xyz, 2));
for patch = 1:npatches
    cols = (patch - 1) * nodes_per_patch + (1:nodes_per_patch);
    R = xyz(:, cols);
    coef = R * A.';
    du = coef * Du;
    dv = coef * Dv;
    normal = cross(du, dv, 1);
    normal_norm = vecnorm(normal, 2, 1);
    bad = normal_norm == 0 | ~isfinite(normal_norm);
    normal(:, ~bad) = normal(:, ~bad) ./ normal_norm(~bad);
    normal(:, bad) = NaN;

    srcvals(1:3, cols) = R;
    srcvals(4:6, cols) = du;
    srcvals(7:9, cols) = dv;
    srcvals(10:12, cols) = normal;
end

norders = order * ones(npatches, 1);
iptype = ones(npatches, 1);
if isfield(mesh, 'iptype')
    iptype(:) = double(mesh.iptype);
end

if isfield(mesh, 'meta')
    meta = mesh.meta;
else
    meta = struct();
end
end

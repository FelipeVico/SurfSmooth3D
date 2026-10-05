function mapping = reference_map(sourceLevel, outputLevel, order, uv)
%REFERENCE_MAP Locate output nodes in the fixed STEP lattice, not in space.
if nargin < 4
    uv = surfsmooth3d.internal.koorn.rv_nodes(order);
end
sourceMaps = surfsmooth3d.edgepreserve.triangle_maps(sourceLevel, 'step');
outputMaps = surfsmooth3d.edgepreserve.triangle_maps(outputLevel, 'smoother');
n = size(uv,2);
coarse = zeros(2, n*size(outputMaps,3));
for k = 1:size(outputMaps,3)
    coarse(:,(k-1)*n+(1:n)) = outputMaps(:,1,k) + ...
        outputMaps(:,2:3,k)*uv;
end
ids = zeros(1,size(coarse,2));
local = zeros(size(coarse));
% Ascending source indices make shared-boundary ownership deterministic.
for k = 1:size(sourceMaps,3)
    remaining = find(ids == 0);
    if isempty(remaining), break; end
    candidate = sourceMaps(:,2:3,k) \ ...
        (coarse(:,remaining)-sourceMaps(:,1,k));
    inside = all(candidate >= -1e-12,1) & sum(candidate,1) <= 1+1e-12;
    hit = remaining(inside);
    ids(hit) = k;
    local(:,hit) = candidate(:,inside);
end
if any(ids == 0)
    error('edgepreserve:unmappedNode', ...
        '%d output nodes have no source reference triangle.', nnz(ids == 0));
end
mapping = struct('coarse_uv',coarse, 'source_patch',ids, ...
    'source_uv',local, 'source_maps',sourceMaps, 'output_maps',outputMaps);
end

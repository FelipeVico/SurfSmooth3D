function Sref = resample_cad(pkg, S, info)
%RESAMPLE_CAD Match fixed CAD polynomials to explicit mixed-depth leaf maps.
if numel(info.parent_ids) ~= S.npatches || size(info.maps,3) ~= S.npatches || ...
        any(S.iptype ~= 1) || any(S.norders ~= S.norders(1))
    error('adaptivesmoother:layout','Invalid adaptive leaf layout.');
end
sourceMaps = surfsmooth3d.edgepreserve.triangle_maps(pkg.source_level,'step');
rv = surfsmooth3d.internal.koorn.rv_nodes(double(S.norders(1)));
n = size(rv,2);
ns = size(sourceMaps,3);
r = zeros(3,S.npts);
for leaf = 1:S.npatches
    parent = info.parent_ids(leaf);
    if parent < 1 || parent > pkg.coarse_triangle_count || parent ~= fix(parent)
        error('adaptivesmoother:parent','Invalid coarse parent ID.');
    end
    map = info.maps(:,:,leaf);
    coarse = map(:,1)+map(:,2:3)*rv;
    assigned = false(1,n);
    for child = 1:ns
        remaining = find(~assigned);
        if isempty(remaining), break; end
        uv = sourceMaps(:,2:3,child) \ (coarse(:,remaining)-sourceMaps(:,1,child));
        inside = all(uv >= -1e-12,1) & sum(uv,1) <= 1+1e-12;
        hit = remaining(inside);
        if isempty(hit), continue; end
        sourcePatch = (parent-1)*ns+child;
        r(:,(leaf-1)*n+hit) = pkg.surface.srccoefs{sourcePatch}(1:3,:)* ...
            surfsmooth3d.internal.koorn.pols(pkg.source_order,uv(:,inside));
        assigned(hit) = true;
    end
    if ~all(assigned)
        error('adaptivesmoother:unmappedNode','Leaf %d has unmapped reference nodes.',leaf);
    end
end
values = surfsmooth3d.edgepreserve.srcvals_from_positions(S,r);
Sref = surfsmooth3d.surfer(S.npatches,S.norders,values,S.iptype);
end

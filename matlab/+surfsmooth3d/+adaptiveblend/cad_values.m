function [r,du,dv] = cad_values(pkg,info,ids,uv)
%CAD_VALUES Evaluate source charts at explicit leaf coordinates, one-sided.
maps = surfsmooth3d.edgepreserve.triangle_maps(pkg.source_level,'step');
ns = size(maps,3); nt = size(uv,2);
r = nan(3,nt); du = r; dv = r;
for leaf = unique(ids(:)).'
    ix = find(ids==leaf); A = info.maps(:,:,leaf);
    coarse = A(:,1)+A(:,2:3)*uv(:,ix);
    assigned = false(1,numel(ix));
    for k = 1:ns
        remaining = find(~assigned);
        if isempty(remaining), break; end
        local = maps(:,2:3,k)\(coarse(:,remaining)-maps(:,1,k));
        inside = all(local>=-1e-12,1) & sum(local,1)<=1+1e-12;
        hit = remaining(inside);
        if isempty(hit), continue; end
        patch = (info.parent_ids(leaf)-1)*ns+k;
        c = pkg.surface.srccoefs{patch}(1:3,:);
        [pol,pu,pv] = surfsmooth3d.internal.koorn.ders(pkg.source_order,local(:,inside));
        J = maps(:,2:3,k)\A(:,2:3);
        r(:,ix(hit)) = c*pol;
        du(:,ix(hit)) = c*(pu*J(1,1)+pv*J(2,1));
        dv(:,ix(hit)) = c*(pu*J(1,2)+pv*J(2,2));
        assigned(hit) = true;
    end
    if ~all(assigned), error('adaptiveblend:unmapped','Unmapped CAD coordinates on leaf %d.',leaf); end
end
end

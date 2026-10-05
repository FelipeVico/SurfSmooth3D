function records = assess_checks(handle,pkg,info,records,ids,b,s,diameter,radius,root)
%ASSESS_CHECKS Independent blended positions and normals, never Phi(B)-1/2.
capacity = max(1,floor(30000/b.nc));
for first = 1:capacity:numel(ids)
    batch = ids(first:min(numel(ids),first+capacity-1));
    owners = repelem(batch(:).',b.nc); uv = repmat(b.check,1,numel(batch));
    [cad,cu,cv] = surfsmooth3d.adaptiveblend.cad_values(pkg,info,owners,uv);
    if isempty(s.selectedEdgeIds)
        reference = cad; refu = cu; refv = cv;
    else
        [smooth,ier] = surfsmooth3d_adaptive_blend_routs('project',handle,owners,uv);
        surfsmooth3d.adaptiveblend.check_error(ier,root,'Independent blended-reference Newton');
        [sg,ier] = surfsmooth3d_adaptive_blend_routs('sigma',handle,cad);
        surfsmooth3d.adaptiveblend.check_error(ier,root,'Independent sigma evaluation');
        [beta,gb] = surfsmooth3d.adaptiveblend.beta_values(cad,pkg.edges,s,sg(1,:),sg(2:4,:),diameter);
        reference = cad+beta.*(smooth(1:3,:)-cad);
        refu = (1-beta).*cu+beta.*smooth(4:6,:)+sum(gb.*cu,1).*(smooth(1:3,:)-cad);
        refv = (1-beta).*cv+beta.*smooth(7:9,:)+sum(gb.*cv,1).*(smooth(1:3,:)-cad);
    end
    for j = 1:numel(batch)
        k = batch(j); ix = (j-1)*b.nc+(1:b.nc); rec = records{k};
        records{k}.indicators(3:4) = surfsmooth3d.adaptiveblend.independent_indicators(rec,b, ...
            reference(:,ix),refu(:,ix),refv(:,ix),info.launch_diameter(k),radius);
        records{k}.checked = true;
    end
end
end

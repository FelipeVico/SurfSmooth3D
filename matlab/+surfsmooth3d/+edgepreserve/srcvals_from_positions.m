function srcvals = srcvals_from_positions(S, rvals)
%SRCVALS_FROM_POSITIONS Differentiate the sampled patch polynomial spectrally.
if any(S.iptype(:) ~= 1) || ~isequal(size(rvals),size(S.r))
    error('edgepreserve:positionLayout', 'Expected matching triangular RV data.');
end
srcvals = zeros(12,S.npts);
srcvals(1:3,:) = rvals;
cache = struct();
for k = 1:S.npatches
    order = double(S.norders(k));
    key = sprintf('p%d',order);
    if ~isfield(cache,key)
        rv = surfsmooth3d.internal.koorn.rv_nodes(order);
        A = surfsmooth3d.internal.koorn.vals2coefs(order,rv);
        [~,Du,Dv] = surfsmooth3d.internal.koorn.ders(order,rv);
        cache.(key) = {A,Du,Dv};
    end
    operators = cache.(key);
    indices = double(S.ixyzs(k)):double(S.ixyzs(k+1))-1;
    coefficients = rvals(:,indices)*operators{1}.';
    du = coefficients*operators{2};
    dv = coefficients*operators{3};
    normal = cross(du,dv,1);
    jac = sqrt(sum(normal.^2,1));
    if any(~isfinite(jac) | jac <= eps)
        error('edgepreserve:degenerateBlendPatch', ...
            'Patch %d has degenerate spectral derivatives.',k);
    end
    srcvals(4:6,indices) = du;
    srcvals(7:9,indices) = dv;
    srcvals(10:12,indices) = normal./jac;
end
end

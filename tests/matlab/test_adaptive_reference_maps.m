function test_adaptive_reference_maps()
%TEST_ADAPTIVE_REFERENCE_MAPS Mixed depths and both source subdivision orders.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
for sourceLevel = 0:2
    sourceMaps = surfsmooth3d.edgepreserve.triangle_maps(sourceLevel,'step');
    sourceIds = repelem((1:2).',size(sourceMaps,3));
    sourceMaps = repmat(sourceMaps,1,1,2);
    source = polynomial_surface(sourceMaps,sourceIds,4);
    pkg = struct('surface',source,'source_level',sourceLevel,'source_order',4,'coarse_triangle_count',2);
    for order = [2 4 6]
        for level = 0:3
            maps = surfsmooth3d.edgepreserve.triangle_maps(level,'smoother');
            % One coarse patch, plus a differently refined second parent.
            maps = cat(3,[0 1 0;0 0 1],maps);
            ids = [1;2*ones(size(maps,3)-1,1)];
            expected = polynomial_surface(maps,ids,order);
            info = struct('maps',maps,'parent_ids',ids);
            sampled = surfsmooth3d.adaptivesmoother.resample_cad(pkg,expected,info);
            assert(max(abs(sampled.r-expected.r),[],'all')<1e-11);
            assert(max(abs(sampled.du-expected.du),[],'all')<1e-11);
            assert(max(abs(sampled.dv-expected.dv),[],'all')<1e-11);
            report = surfsmooth3d.edgepreserve.validity_report(sampled,'polynomial',expected);
            assert(report.is_valid && max(abs(report.jac_ratio_to_ref-1))<1e-10);
        end
    end
end
fprintf('PASS: mixed-depth CAD correspondence, source levels 0:2, output depths 0:3 and derivative scaling.\n');
end

function S = polynomial_surface(maps,ids,order)
rv = surfsmooth3d.internal.koorn.rv_nodes(order); n = size(rv,2); nt = numel(ids);
r = zeros(3,n*nt);
for k=1:nt
    uv = maps(:,1,k)+maps(:,2:3,k)*rv;
    u=uv(1,:); v=uv(2,:);
    r(:,(k-1)*n+(1:n)) = [u+10*ids(k);v+ids(k);ids(k)+.1*u.^2+.07*u.*v+.08*v.^2];
end
template = struct('npatches',nt,'npts',n*nt,'norders',order*ones(nt,1), ...
    'iptype',ones(nt,1),'r',r,'ixyzs',(1:n:n*nt+1).');
S = surfsmooth3d.surfer(nt,order,surfsmooth3d.edgepreserve.srcvals_from_positions(template,r));
end

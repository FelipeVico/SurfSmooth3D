function test_adaptive_blend_math()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
s = surfsmooth3d.adaptiveblend.default_settings(); s.selectedEdgeIds = 1;
x = linspace(-2,2,101); w = ones(size(x))*.04; w([1 end]) = .02;
edges = struct('edge_id',1,'panels',struct('xyz',[x;0*x;0*x],'w',w));
r = [.1 .5 -.7;.2 .4 .7;.3 .1 -.2];
g = [.013;-.007;.009]; sigma = .2+g.'*r;
[beta,gb] = surfsmooth3d.adaptiveblend.beta_values(r,edges,s,sigma,repmat(g,1,3),4);
for axis = 1:3
    dr = zeros(size(r)); dr(axis,:) = 1e-6;
    plus = surfsmooth3d.adaptiveblend.beta_values(r+dr,edges,s,.2+g.'*(r+dr),repmat(g,1,3),4);
    minus = surfsmooth3d.adaptiveblend.beta_values(r-dr,edges,s,.2+g.'*(r-dr),repmat(g,1,3),4);
    assert(max(abs((plus-minus)/2e-6-gb(axis,:)))<2e-8);
end
assert(all(beta>0 & beta<1));
% Empty selections and saturated/floored distances retain the old formula.
empty = s; empty.selectedEdgeIds = [];
[be,ge] = surfsmooth3d.adaptiveblend.beta_values(r,edges,empty,sigma,repmat(g,1,3),4);
assert(all(be==0) && all(ge==0,'all'));
for p = 1:6
    b = surfsmooth3d.adaptiveblend.basis(p);
    assert(b.nc==4*b.np+3*(p+2)+3);
    r = [b.rv;zeros(1,b.np)];
    rec = surfsmooth3d.adaptiveblend.nodal_patch(r,r,zeros(1,b.np),ones(1,b.np),b,sqrt(2),1);
    assert(max(rec.indicators(1:2))<1e-12);
    rec = surfsmooth3d.adaptiveblend.nodal_patch(r*0,r,zeros(1,b.np),ones(1,b.np),b,sqrt(2),1);
    assert(all(isinf(rec.indicators(1:2))));
end
% A narrow bump at an independent interior point is invisible at the RV grid.
b = surfsmooth3d.adaptiveblend.basis(4); center = b.check(:,2); width = .003;
bump = @(uv) .03*exp(-sum((uv-center).^2,1)/width^2);
r = [b.rv;bump(b.rv)];
rec = surfsmooth3d.adaptiveblend.nodal_patch(r,r,zeros(1,b.np),ones(1,b.np),b,sqrt(2),1);
q = bump(b.check); ref = [b.check;q];
ru = [ones(1,b.nc);zeros(1,b.nc);-2*(b.check(1,:)-center(1)).*q/width^2];
rv = [zeros(1,b.nc);ones(1,b.nc);-2*(b.check(2,:)-center(2)).*q/width^2];
eta = surfsmooth3d.adaptiveblend.independent_indicators(rec,b,ref,ru,rv,sqrt(2),1);
assert(max(rec.indicators(1:2))<1e-4 && eta(1)>.01);
etaFlip = surfsmooth3d.adaptiveblend.independent_indicators(rec,b,ref,rv,ru,sqrt(2),1);
assert(all(isinf(etaFlip)));
% Both comparisons are invariant under translation and uniform scaling.
scale = 13; shift = [17;-4;8];
rec2 = surfsmooth3d.adaptiveblend.nodal_patch(scale*r+shift,scale*r+shift, ...
    zeros(1,b.np),ones(1,b.np),b,scale*sqrt(2),scale);
eta2 = surfsmooth3d.adaptiveblend.independent_indicators(rec2,b,scale*ref+shift,scale*ru,scale*rv,scale*sqrt(2),scale);
assert(max(abs(eta-eta2))<1e-11);
test_cad_maps();
files = dir(fullfile(root,'matlab','+surfsmooth3d','+adaptiveblend','*.m'));
for k = 1:numel(files)
    messages = checkcode(fullfile(files(k).folder,files(k).name),'-id');
    for j = 1:numel(messages)
        fprintf('%s:%d %s %s\n',files(k).name,messages(j).line,messages(j).id,messages(j).message);
    end
    assert(~any(strcmp({messages.id},'PARSE')),'MATLAB parse error.');
end
fprintf('PASS: adaptive blend beta gradients, planar tails, low orders and degeneracy.\n');
end

function test_cad_maps()
for sourceLevel = 0:2
    maps = surfsmooth3d.edgepreserve.triangle_maps(sourceLevel,'step');
    rv = surfsmooth3d.internal.koorn.rv_nodes(4); np = size(rv,2); n = size(maps,3);
    r = zeros(3,np*n);
    for k = 1:n
        uv = maps(:,1,k)+maps(:,2:3,k)*rv;
        r(:,(k-1)*np+(1:np)) = poly(uv);
    end
    template = struct('npatches',n,'npts',n*np,'norders',4*ones(n,1), ...
        'iptype',ones(n,1),'r',r,'ixyzs',(1:np:n*np+1).');
    source = surfsmooth3d.surfer(n,4,surfsmooth3d.edgepreserve.srcvals_from_positions(template,r),1);
    pkg = struct('surface',source,'source_level',sourceLevel,'source_order',4);
    for level = 0:3
        maps = surfsmooth3d.edgepreserve.triangle_maps(level,'smoother');
        info = struct('maps',maps,'parent_ids',ones(size(maps,3),1));
        uv = [0 1 0 .5 .5 0 .21;0 0 1 0 .5 .5 .32];
        for k = 1:size(maps,3)
            [r,du,dv] = surfsmooth3d.adaptiveblend.cad_values(pkg,info,k*ones(1,size(uv,2)),uv);
            coarse = maps(:,1,k)+maps(:,2:3,k)*uv;
            [expected,cu,cv] = poly(coarse);
            J = maps(:,2:3,k);
            assert(max(abs(r-expected),[],'all')<1e-11);
            assert(max(abs(du-cu*J(1,1)-cv*J(2,1)),[],'all')<1e-11);
            assert(max(abs(dv-cu*J(1,2)-cv*J(2,2)),[],'all')<1e-11);
        end
    end
end
end

function [r,du,dv] = poly(uv)
u = uv(1,:); v = uv(2,:);
r = [u+4;v-2;.1*u.^2+.07*u.*v+.08*v.^2];
du = [ones(size(u));zeros(size(u));.2*u+.07*v];
dv = [zeros(size(u));ones(size(u));.07*u+.16*v];
end

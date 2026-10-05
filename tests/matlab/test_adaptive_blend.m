function test_adaptive_blend(mode)
%TEST_ADAPTIVE_BLEND Native sessions, compatibility, selective solves and export.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder,'s'));
pkg = sphere_package(folder);
lockBaseline = adaptive_blend_lock_baseline();
if nargin>0 && strcmp(mode,'limits')
    test_limits(pkg,lockBaseline);
    return
end
s = surfsmooth3d.adaptiveblend.default_settings(); s.order = 4; s.max_refine = 0; s.eps_adapt = .01;
cadfile = fullfile(folder,'source.txt'); surfsmooth3d.edgepreserve.write_point_source(pkg.surface,cadfile);
opts = struct('fcad',cadfile,'nquad',8,'filetype',3,'rlam',2, ...
    'two_stage_smoother',true,'adapt_sigma',1,'nrefine',1);
uniform = surfsmooth3d.multiscale_mesher(pkg.scaffold_file,4,opts);
outroot = fullfile(folder,'native');
[h,ier] = surfsmooth3d_adaptive_blend_routs('open',pkg.scaffold_file,cadfile,outroot,8,4,1,2,200000);
assert(ier==0); session = onCleanup(@() surfsmooth3d_adaptive_blend_routs('close',h));
[values,meta,sphere] = surfsmooth3d_adaptive_blend_routs('get',h);
assert(max(abs(values-extract_arrays(uniform{1})),[],'all')<2e-11);
b = surfsmooth3d.adaptiveblend.basis(4);
owners = repelem([1 7],b.np); uv = repmat(b.rv,1,2);
[projected,ier] = surfsmooth3d_adaptive_blend_routs('project',h,owners,uv); assert(ier==0);
nodes = [1:b.np,6*b.np+(1:b.np)];
assert(max(abs(projected-values(:,nodes)),[],'all')<2e-8);
[sigma0,ier] = surfsmooth3d_adaptive_blend_routs('sigma',h,pkg.surface.r(:,1:91:end)); assert(ier==0);
assert(sphere(4)>0 && size(meta,2)==pkg.surface.npatches);
ier = surfsmooth3d_adaptive_blend_routs('refine',h,[1 7],size(values,2)); assert(ier==8);
[unchanged,~,~] = surfsmooth3d_adaptive_blend_routs('get',h); assert(isequal(values,unchanged));
ier = surfsmooth3d_adaptive_blend_routs('refine',h,[1 7],200000); assert(ier==0);
[children,maps,~] = surfsmooth3d_adaptive_blend_routs('get',h);
assert(size(maps,2)==size(meta,2)+6);
assert(isequal(children(:,4*b.np+(1:b.np)),values(:,b.np+(1:b.np))));
childMap = surfsmooth3d.edgepreserve.triangle_maps(1,'smoother');
assert(max(abs(reshape(maps(3:8,1:4),2,3,4)-childMap),[],'all')<1e-14);
u = extract_arrays(uniform{2});
assert(max(abs(children(:,1:4*b.np)-u(:,1:4*b.np)),[],'all')<2e-8);
[sigma1,ier] = surfsmooth3d_adaptive_blend_routs('sigma',h,pkg.surface.r(:,1:91:end));
assert(ier==0 && isequal(sigma0,sigma1));
clear session
must_fail(@() stale(h));
empty = surfsmooth3d.adaptiveblend.solve(pkg,s);
assert(all(empty.beta==0));
assert(max(abs(empty.surface.r-empty.reference.r),[],'all')==0);
s.selectedEdgeIds = 1;
result = surfsmooth3d.adaptiveblend.solve(pkg,s);
params = s; params.sigma = result.sigma;
params.geometryDiameter = norm(max(pkg.surface.r,[],2)-min(pkg.surface.r,[],2));
[old,beta] = surfsmooth3d.edgepreserve.blend(result.reference,result.smooth,pkg.edges,params);
assert(max(abs(beta-result.beta))<2e-14);
assert(max(abs(extract_arrays(old)-extract_arrays(result.surface)),[],'all')<2e-9);
assert(isequal(empty.smooth.r,result.smooth.r));
s.max_refine = 2; s.eps_adapt = 3e-3;
refined = surfsmooth3d.adaptiveblend.solve(pkg,s);
assert(refined.surface.npatches>result.surface.npatches);
assert(refined.surface.npatches<result.surface.npatches*16);
assert(any(refined.info.depth<refined.info.achieved_depth));
assert(refined.report.is_valid);
manifest = surfsmooth3d.adaptiveblend.save_result(pkg,refined,folder);
assert(all(manifest.status=="saved"));
loaded = surfsmooth3d.surfer.load_from_file(char(manifest.filename(1)));
assert(max(abs(extract_arrays(loaded)-extract_arrays(refined.surface)),[],'all')<1e-11);
parameters = load(fullfile(fileparts(manifest.filename(1)),'parameters.mat'));
assert(isequal(parameters.parameters.info.maps,refined.info.maps));
refined.report.is_valid = false;
manifest = surfsmooth3d.adaptiveblend.save_result(pkg,refined,folder);
assert(manifest.status(1)=="rejected" && manifest.status(2)=="saved");
refined.report_smooth.is_valid = false;
must_fail(@() surfsmooth3d.adaptiveblend.save_result(pkg,refined,folder));
assert(mislocked('surfsmooth3d_adaptive_blend_routs')==lockBaseline);
for k = 1:numel(pkg.source_paths)
    assert(strcmp(pkg.source_hashes{k},surfsmooth3d.edgepreserve.file_sha256(pkg.source_paths{k})));
end
fprintf('PASS: adaptive blend native/session, baseline, selective refinement, fixed sigma and verified pair export.\n');
end

function test_limits(pkg,lockBaseline)
s = surfsmooth3d.adaptiveblend.default_settings(); s.order = 8; s.max_refine = 0; s.eps_adapt = .01;
result = surfsmooth3d.adaptiveblend.solve(pkg,s);
assert(result.info.converged && all(result.info.independent_checked));
s.eps_adapt = 1e-12; s.max_refine = 2; s.max_points = result.surface.npts;
limited = surfsmooth3d.adaptiveblend.solve(pkg,s);
assert(strcmp(limited.info.stop_reason,'point budget reached'));
assert(isequal(limited.surface.r,result.surface.r));
s.max_points = s.max_points-1;
must_fail(@() surfsmooth3d.adaptiveblend.solve(pkg,s));
assert(mislocked('surfsmooth3d_adaptive_blend_routs')==lockBaseline);
s = surfsmooth3d.adaptiveblend.default_settings(); s.order = 4; s.max_refine = 2;
s.eps_adapt = .01; s.selectedEdgeIds = 1; s.c_rad = 100;
full = surfsmooth3d.adaptiveblend.solve(pkg,s);
assert(all(full.beta==1) && full.info.converged && full.report.is_valid);
assert(max(abs(full.surface.r-full.smooth.r),[],'all')==0);
fprintf('PASS: tolerance achieved, order-eight CAD-only case, forced beta=1 and point-limit handling.\n');
end

function stale(h)
[~,~,~] = surfsmooth3d_adaptive_blend_routs('get',h);
end
function must_fail(call)
try
    call();
catch
    return
end
error('adaptiveblend:test','Expected failure.');
end
function pkg = sphere_package(folder)
S = fixture_sphere(1,2,[0;0;0],8,1);
pol = surfsmooth3d.internal.koorn.pols(8,[0 1 0;0 0 1]); points = zeros(3,3*S.npatches);
for k = 1:S.npatches, points(:,3*k-2:3*k) = S.srccoefs{k}(1:3,:)*pol; end
[nodes,~,idx] = unique(round(points.',10),'rows'); tris = reshape(idx,3,[]);
fname = fullfile(folder,'scaffold.gidmsh');
fid = fopen(fname,'w');
fprintf(fid,'MESH dimension 3 ElemType Triangle Nnode 3\nCoordinates\n');
fprintf(fid,'%d %.17g %.17g %.17g\n',[(1:size(nodes,1)).' nodes].');
fprintf(fid,'End Coordinates\nElements\n');
fprintf(fid,'%d %d %d %d 1\n',[(1:size(tris,2));tris]); fprintf(fid,'End Elements\n'); fclose(fid);
t = (0:199)*2*pi/200;
edge = struct('edge_id',1,'panels',struct('xyz',[cos(t);sin(t);0*t],'w',ones(1,200)*2*pi/200));
pkg = struct('surface',S,'source_order',8,'source_level',0,'coarse_triangle_count',S.npatches, ...
    'scaffold_file',fname,'source_paths',{{fname}},'source_hashes',{{surfsmooth3d.edgepreserve.file_sha256(fname)}}, ...
    'edges',edge,'name','sphere','geometry_name','sphere');
end

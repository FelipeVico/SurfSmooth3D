function test_adaptive_smoother()
%TEST_ADAPTIVE_SMOOTHER End-to-end baseline, fixed sigma, limits and export.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
directory = tempname; mkdir(directory);
cleanup = onCleanup(@() rmdir(directory,'s'));
source = fixture_sphere(1,2,[0;0;0],8,1);
scaffold = fullfile(directory,'scaffold.gidmsh');
write_scaffold(source,scaffold);
cad = fullfile(directory,'cad_points.txt');
surfsmooth3d.edgepreserve.write_point_source(source,cad);
hash = surfsmooth3d.edgepreserve.file_sha256(scaffold);
opts = struct('fcad',cad,'filetype',3,'nquad',8,'rlam',2, ...
    'adapt_sigma',1,'two_stage_smoother',true,'nrefine',2, ...
    'max_refine',0,'eps_adapt',1e-12,'max_points',200000);
uniform = surfsmooth3d.multiscale_mesher(scaffold,4,opts);
[baseline,zero] = surfsmooth3d.multiscale_mesher_adaptive(scaffold,4,opts);
assert_same(baseline,uniform{1},2e-11);
assert(zero.achieved_depth==0 && ~zero.converged);
opts.max_refine = 2;
[adaptive,info] = surfsmooth3d.multiscale_mesher_adaptive(scaffold,4,opts);
assert(all(info.depth==2),'Forced-all test did not refine every leaf.');
assert_same(adaptive,uniform{3},2e-9);
assert(adaptive.npatches==source.npatches*16);
maps = surfsmooth3d.edgepreserve.triangle_maps(2,'smoother');
assert(max(abs(info.maps(:,:,1:16)-maps),[],'all')<1e-14);
assert(all(info.parent_ids==repelem((1:source.npatches).',16)));
determinants = squeeze(info.maps(1,2,:).*info.maps(2,3,:)-info.maps(2,2,:).*info.maps(1,3,:));
assert(all(abs(determinants-4^-2)<1e-14));
opts.max_points = baseline.npts;
[limited,limitInfo] = surfsmooth3d.multiscale_mesher_adaptive(scaffold,4,opts);
assert_same(limited,baseline,2e-11);
assert(strcmp(limitInfo.stop_reason,'point budget reached'));
opts.max_points = baseline.npts-1;
must_fail(@() surfsmooth3d.multiscale_mesher_adaptive(scaffold,4,opts));
opts.max_points = 200000;
opts.max_refine = 0;
for mode = [0 1 3]
    opts.adapt_sigma = mode;
    targets = source.r(:,1:71:end);
    [sigma0,gradient0] = surfsmooth3d.multiscale_mesher_sigma_eval(scaffold,4,opts,targets);
    surfsmooth3d.multiscale_mesher_adaptive(scaffold,4,opts);
    opts.nrefine = 3;
    [sigma1,gradient1] = surfsmooth3d.multiscale_mesher_sigma_eval(scaffold,8,opts,targets);
    assert(isequal(sigma0,sigma1) && isequal(gradient0,gradient1));
end
assert(strcmp(hash,surfsmooth3d.edgepreserve.file_sha256(scaffold)));
opts.eps_adapt = .99;
[lowOrder,lowInfo] = surfsmooth3d.multiscale_mesher_adaptive(scaffold,1,opts);
assert(all(lowOrder.norders==1) && all(lowInfo.independent_checked));
assert(all(lowInfo.indicators(:,1)==0),'Order-one position tails must be empty.');
% Mixed-level maps use independent source patches, never physical proximity.
pkg = struct('surface',source,'source_level',0,'source_order',8, ...
    'coarse_triangle_count',source.npatches);
reference = surfsmooth3d.adaptivesmoother.resample_cad(pkg,adaptive,info);
report = surfsmooth3d.edgepreserve.validity_report(adaptive,'test',reference);
assert(report.is_valid);
pkg.name = 'sphere'; pkg.geometry_name = 'sphere';
pkg.source_paths = {scaffold,cad}; pkg.source_hashes = {hash,surfsmooth3d.edgepreserve.file_sha256(cad)};
settings = struct('order',4,'rlam',2,'adapt_sigma',1,'eps_adapt',1e-12, ...
    'max_refine',2,'max_points',200000);
result = struct('surface',adaptive,'info',info,'report',report,'settings',settings);
filename = surfsmooth3d.adaptivesmoother.save_result(pkg,result,directory);
assert(contains(filename,'_tolunmet') && isfile(filename));
assert_same(surfsmooth3d.surfer.load_from_file(filename),adaptive,1e-12);
saved = load(fullfile(fileparts(filename),'parameters.mat'));
assert(isequal(saved.parameters.info.maps,info.maps));
result.report.is_valid = false;
must_fail(@() surfsmooth3d.adaptivesmoother.save_result(pkg,result,directory));
fprintf('PASS: adaptive baseline/uniform equivalence, maps, fixed sigma, limits and verified export.\n');
end

function assert_same(a,b,tolerance)
assert(a.npatches==b.npatches && isequal(a.norders,b.norders));
x = extract_arrays(a); y = extract_arrays(b);
assert(max(abs(x(:)-y(:)))<tolerance*max(1,max(abs(y(:)))), ...
    'Adaptive/uniform data mismatch: %.4g.',max(abs(x(:)-y(:))));
end

function must_fail(call)
try
    call();
catch
    return
end
error('adaptive_test:expectedFailure','Expected operation to fail.');
end

function write_scaffold(S,filename)
pol = surfsmooth3d.internal.koorn.pols(double(S.norders(1)),[0 1 0;0 0 1]);
points = zeros(3,3*S.npatches);
for k=1:S.npatches
    points(:,3*k-2:3*k) = S.srccoefs{k}(1:3,:)*pol;
end
[nodes,~,indices] = unique(round(points.',10),'rows');
tris = reshape(indices,3,[]);
fid = fopen(filename,'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'MESH dimension 3 ElemType Triangle Nnode 3\nCoordinates\n');
fprintf(fid,'%d %.17g %.17g %.17g\n',[(1:size(nodes,1)).' nodes].');
fprintf(fid,'End Coordinates\nElements\n');
fprintf(fid,'%d %d %d %d 1\n',[(1:size(tris,2));tris]);
fprintf(fid,'End Elements\n');
end

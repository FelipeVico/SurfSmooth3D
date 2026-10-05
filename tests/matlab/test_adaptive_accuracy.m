function test_adaptive_accuracy(savedResult,packageDir)
%TEST_ADAPTIVE_ACCURACY Compare mixed-depth charts with a finer fixed-source solve.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
if nargin < 1, savedResult = fullfile(tempdir,'adaptive_gui_final/adaptive_real_result.mat'); end
if nargin < 2
    packageDir = fullfile(fullfile(root,'tests','fixtures','cad'), ...
        'convergence_intersection_two_balls','p4_mf0p10000_rl0');
end
stored = load(savedResult); tight = stored.result;
pkg = surfsmooth3d.edgepreserve.load_package(packageDir);
[work,cleanup] = surfsmooth3d.edgepreserve.stage_source(pkg); %#ok<ASGLU>
opts = struct('nquad',pkg.source_order,'fcad',work.cad_file,'filetype',3, ...
    'rlam',tight.settings.rlam,'adapt_sigma',tight.settings.adapt_sigma, ...
    'two_stage_smoother',true,'nrefine',2);
reference = surfsmooth3d.multiscale_mesher(work.scaffold_file,6,opts);
reference = reference{end};
looseSettings = tight.settings; looseSettings.eps_adapt = 3e-3;
loose = surfsmooth3d.adaptivesmoother.solve(pkg,looseSettings);
errors = zeros(2,3);
cases = {loose,tight};
for k=1:2
    errors(k,1:2) = compare_reference(cases{k},reference,2);
    errors(k,3) = shared_edges(cases{k},pkg);
end
fprintf('Independent position / normal / shared-edge errors (loose, then tight):\n');
fprintf('  %.6e  %.6e  %.6e\n',errors.');
assert(all(errors(2,1:2)<errors(1,1:2)),'Tighter tolerance did not improve measured geometric error.');
assert(errors(2,1)<1e-3 && errors(2,2)<3e-3 && errors(2,3)<1e-3);
assert(tight.surface.npatches < pkg.coarse_triangle_count*4^tight.info.achieved_depth);
fprintf('PASS: mixed-depth accuracy improves with tolerance; shared-edge discrepancies measured.\n');
end

function errors = compare_reference(result,reference,level)
S = result.surface; info = result.info;
uv = [surfsmooth3d.internal.koorn.rv_nodes(7) [0 .25 .5 .75 1 0;0 0 0 0 0 1]];
[pol,du,dv] = surfsmooth3d.internal.koorn.ders(double(S.norders(1)),uv);
maps = surfsmooth3d.edgepreserve.triangle_maps(level,'smoother');
errors = [0 0];
for k=1:S.npatches
    coarse = info.maps(:,1,k)+info.maps(:,2:3,k)*uv;
    r = S.srccoefs{k}(1:3,:)*pol;
    n = cross(S.srccoefs{k}(1:3,:)*du,S.srccoefs{k}(1:3,:)*dv,1);
    n = n./vecnorm(n);
    remaining = true(1,size(uv,2));
    for child=1:size(maps,3)
        indices = find(remaining);
        local = maps(:,2:3,child) \ (coarse(:,indices)-maps(:,1,child));
        inside = all(local>=-1e-12,1) & sum(local,1)<=1+1e-12;
        hit = indices(inside);
        if isempty(hit), continue; end
        patch = (info.parent_ids(k)-1)*size(maps,3)+child;
        [rp,ru,rv] = surfsmooth3d.internal.koorn.ders(double(reference.norders(1)),local(:,inside));
        coefficients = reference.srccoefs{patch}(1:3,:);
        expected = coefficients*rp;
        normal = cross(coefficients*ru,coefficients*rv,1);
        normal = normal./vecnorm(normal);
        errors = max(errors,[max(vecnorm(r(:,hit)-expected))/info.radius,max(vecnorm(n(:,hit)-normal))]);
        remaining(hit) = false;
    end
    assert(~any(remaining));
end
end

function discrepancy = shared_edges(result,pkg)
S = result.surface; info = result.info;
t=linspace(0,1,17); uv=[t 1-t zeros(size(t));zeros(size(t)) t 1-t];
pol=surfsmooth3d.internal.koorn.pols(double(S.norders(1)),uv);
launch=zeros(3,size(uv,2)*S.npatches); values=launch;
for k=1:S.npatches
    coarse=info.maps(:,1,k)+info.maps(:,2:3,k)*uv;
    triangle=pkg.scaffold.nodes(:,pkg.scaffold.tris(:,info.parent_ids(k)));
    indices=(k-1)*size(uv,2)+(1:size(uv,2));
    launch(:,indices)=triangle*[1-sum(coarse,1);coarse];
    values(:,indices)=S.srccoefs{k}(1:3,:)*pol;
end
[~,~,group]=unique(round(launch.'/info.radius,11),'rows');
spread=zeros(max(group),3);
for axis=1:3
    spread(:,axis)=accumarray(group,values(axis,:).',[],@max)-accumarray(group,values(axis,:).',[],@min);
end
discrepancy=max(vecnorm(spread,2,2))/info.radius;
end

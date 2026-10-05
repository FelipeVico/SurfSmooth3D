function test_adaptive_blend_accuracy(cachedResult)
%TEST_ADAPTIVE_BLEND_ACCURACY Common independent samples of one fixed blend.
% Optionally reuse result.mat from smoke_adaptive_blend_gui.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
if nargin<1 || isempty(cachedResult)
    packageDir = fullfile(fullfile(root,'tests','fixtures','cad'), ...
        'convergence_intersection_two_balls','p4_mf0p10000_rl0');
    pkg = surfsmooth3d.edgepreserve.load_package(packageDir);
    settings = surfsmooth3d.adaptiveblend.default_settings();
    settings.order = 4; settings.max_refine = 2; settings.eps_adapt = 3e-3;
    settings.selectedEdgeIds = [pkg.edges.edge_id];
    tight = surfsmooth3d.adaptiveblend.solve(pkg,settings);
else
    saved = load(cachedResult); pkg = saved.pkg; tight = saved.result;
end
settings = tight.settings; settings.eps_adapt = .02;
loose = surfsmooth3d.adaptiveblend.solve(pkg,settings);
[work,cleanup] = surfsmooth3d.edgepreserve.stage_source(pkg); %#ok<ASGLU>
outroot = fullfile(work.directory,'accuracy');
[h,ier] = surfsmooth3d_adaptive_blend_routs('open',work.scaffold_file,work.cad_file, ...
    outroot,pkg.source_order,settings.order,settings.adapt_sigma,settings.rlam,settings.max_points);
surfsmooth3d.adaptiveblend.check_error(ier,outroot,'Accuracy reference initialization');
session = onCleanup(@() surfsmooth3d_adaptive_blend_routs('close',h));
% A fixed check grid, unrelated to either adaptive leaf layout.
uv = surfsmooth3d.internal.koorn.rv_nodes(7); np = size(uv,2);
owners = repelem(1:pkg.coarse_triangle_count,np);
uv = repmat(uv,1,pkg.coarse_triangle_count);
info = struct('maps',repmat([0 1 0;0 0 1],1,1,pkg.coarse_triangle_count), ...
    'parent_ids',(1:pkg.coarse_triangle_count).');
[cad,cu,cv] = surfsmooth3d.adaptiveblend.cad_values(pkg,info,owners,uv);
[smooth,ier] = surfsmooth3d_adaptive_blend_routs('project',h,owners,uv);
surfsmooth3d.adaptiveblend.check_error(ier,outroot,'Accuracy reference projection');
[sg,ier] = surfsmooth3d_adaptive_blend_routs('sigma',h,cad);
surfsmooth3d.adaptiveblend.check_error(ier,outroot,'Accuracy sigma');
diameter = max(norm(max(pkg.surface.r,[],2)-min(pkg.surface.r,[],2)),1);
[beta,gb] = surfsmooth3d.adaptiveblend.beta_values(cad,pkg.edges,settings,sg(1,:),sg(2:4,:),diameter);
ref = cad+beta.*(smooth(1:3,:)-cad);
du = (1-beta).*cu+beta.*smooth(4:6,:)+sum(gb.*cu,1).*(smooth(1:3,:)-cad);
dv = (1-beta).*cv+beta.*smooth(7:9,:)+sum(gb.*cv,1).*(smooth(1:3,:)-cad);
normal = cross(du,dv,1); normal = normal./vecnorm(normal);
errors = zeros(2,4);
cases = {loose,tight};
for j = 1:2
    [r,n] = evaluate(cases{j},owners,uv);
    ep = vecnorm(r-ref)/diameter; en = vecnorm(n-normal);
    assert(all(isfinite(ep)) && all(isfinite(en)));
    errors(j,:) = [sqrt(mean(ep.^2)),sqrt(mean(en.^2)),max(ep),max(en)];
end
disp(array2table(errors,'VariableNames',{'position_rms','normal_rms','position_max','normal_max'}, ...
    'RowNames',{'loose','tight'}));
assert(all(errors(2,1:2)<errors(1,1:2)),'Tightening tolerance should reduce independent RMS errors.');
fprintf('PASS: partial-surface position and normal RMS errors decrease on fixed independent samples.\n');
end

function [r,n] = evaluate(result,owners,uv)
r = nan(3,size(uv,2)); n = r;
for coarse = unique(owners)
    indices = find(owners==coarse); assigned = false(size(indices));
    leaves = find(result.info.parent_ids==coarse).';
    for leaf = leaves
        remaining = find(~assigned); if isempty(remaining), break; end
        map = result.info.maps(:,:,leaf);
        local = map(:,2:3)\(uv(:,indices(remaining))-map(:,1));
        inside = all(local>=-1e-12,1) & sum(local,1)<=1+1e-12;
        hit = remaining(inside); if isempty(hit), continue; end
        [pol,pu,pv] = surfsmooth3d.internal.koorn.ders(result.settings.order,local(:,inside));
        c = result.surface.srccoefs{leaf}(1:3,:);
        r(:,indices(hit)) = c*pol;
        normal = cross(c*pu,c*pv,1);
        n(:,indices(hit)) = normal./vecnorm(normal);
        assigned(hit) = true;
    end
    assert(all(assigned));
end
end

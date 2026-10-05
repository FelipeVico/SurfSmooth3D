function result = solve(pkg,s)
%SOLVE Refine the partial surface while keeping the CAD and sigma source fixed.
surfsmooth3d.adaptiveblend.validate_settings(s,pkg);
[work,workCleanup] = surfsmooth3d.edgepreserve.stage_source(pkg); %#ok<ASGLU>
root = fullfile(work.directory,'partial_adaptive');
[handle,ier] = surfsmooth3d_adaptive_blend_routs('open',work.scaffold_file,work.cad_file, ...
    root,pkg.source_order,s.order,s.adapt_sigma,s.rlam,s.max_points);
surfsmooth3d.adaptiveblend.check_error(ier,root,'Initial two-stage solve');
sessionCleanup = onCleanup(@() surfsmooth3d_adaptive_blend_routs('close',handle));
b = surfsmooth3d.adaptiveblend.basis(s.order);
diameter = max(norm(max(pkg.surface.r,[],2)-min(pkg.surface.r,[],2)),1);
records = {}; history = zeros(0,11);
for pass = 0:s.max_refine
    [smooth,meta,sphere] = surfsmooth3d_adaptive_blend_routs('get',handle);
    info = struct('parent_ids',meta(1,:).','depth',meta(2,:).', ...
        'maps',reshape(meta(3:8,:),2,3,[]),'launch_diameter',meta(9,:).');
    n = size(meta,2);
    if isempty(records), records = cell(n,1); end
    fresh = find(cellfun(@isempty,records));
    capacity = max(1,floor(30000/b.np));
    for first = 1:capacity:numel(fresh)
        batch = fresh(first:min(numel(fresh),first+capacity-1));
        owners = repelem(batch(:).',b.np); uv = repmat(b.rv,1,numel(batch));
        cad = surfsmooth3d.adaptiveblend.cad_values(pkg,info,owners,uv);
        [sg,ier] = surfsmooth3d_adaptive_blend_routs('sigma',handle,cad);
        surfsmooth3d.adaptiveblend.check_error(ier,root,'Nodal sigma evaluation');
        beta = surfsmooth3d.adaptiveblend.beta_values(cad,pkg.edges,s,sg(1,:),sg(2:4,:),diameter);
        for j = 1:numel(batch)
            k = batch(j); ix = (j-1)*b.np+(1:b.np); nodes = (k-1)*b.np+(1:b.np);
            r = beta(ix).*smooth(1:3,nodes)+(1-beta(ix)).*cad(:,ix);
            records{k} = surfsmooth3d.adaptiveblend.nodal_patch(r,cad(:,ix),beta(ix),sg(1,ix), ...
                b,info.launch_diameter(k),sphere(4)); %#ok<AGROW> Leaf count is preallocated above.
        end
    end
    eta = cell2mat(cellfun(@(x) x.indicators,records,'UniformOutput',false));
    checked = cellfun(@(x) x.checked,records);
    checkIds = find(~checked & max(eta(:,1:2),[],2)<=s.eps_adapt);
    records = surfsmooth3d.adaptiveblend.assess_checks(handle,pkg,info,records,checkIds,b,s,diameter,sphere(4),root);
    eta = cell2mat(cellfun(@(x) x.indicators,records,'UniformOutput',false));
    checked = cellfun(@(x) x.checked,records);
    unresolved = max(eta,[],2)>s.eps_adapt;
    marked = find(unresolved & info.depth<s.max_refine);
    maxima = max(eta,[],1);
    history(end+1,:) = [pass n n*b.np numel(marked) sum(unresolved) max(info.depth) maxima sum(checked)]; %#ok<AGROW>
    fprintf('Partial adaptive pass %d: %d patches, %d unresolved, depth %d\n', ...
        pass,n,sum(unresolved),max(info.depth));
    fprintf('  max [position tail, normal tail, position check, normal check] = %s\n',mat2str(maxima,4));
    drawnow limitrate nocallbacks;
    reason = 'tolerance achieved';
    if ~any(unresolved), break; end
    reason = 'maximum depth reached';
    if isempty(marked), break; end
    reason = 'point budget reached';
    if n+3*numel(marked)>floor(s.max_points/b.np), break; end
    ier = surfsmooth3d_adaptive_blend_routs('refine',handle,marked,s.max_points);
    surfsmooth3d.adaptiveblend.check_error(ier,root,'Selected-child Newton');
    % Native refinement publishes all four children together, in parent order.
    next = cell(n+3*numel(marked),1); cursor = 0;
    mask = false(n,1); mask(marked) = true;
    for k = 1:n
        if mask(k), cursor = cursor+4; else, cursor = cursor+1; next{cursor} = records{k}; end
    end
    records = next;
end
info.indicators = eta; info.independent_checked = checked;
info.unresolved = unresolved; info.unresolved_count = sum(unresolved);
info.converged = ~any(unresolved); info.stop_reason = reason;
info.achieved_depth = max(info.depth); info.center = sphere(1:3); info.radius = sphere(4);
info.history = array2table(history,'VariableNames',{'pass','patches','nodes','marked', ...
    'unresolved','depth','position_tail','normal_tail','position_check','normal_check','checked'});
info.tolerance_applies_to = 'partial surface only';
Ssmooth = surfsmooth3d.surfer(n,s.order,smooth,1);
src = cellfun(@(x) x.src,records,'UniformOutput',false);
cad = cellfun(@(x) x.cad,records,'UniformOutput',false);
partial = surfsmooth3d.surfer(n,s.order,cat(2,src{:}),1);
reference = surfsmooth3d.surfer(n,s.order,surfsmooth3d.edgepreserve.srcvals_from_positions(Ssmooth,cat(2,cad{:})),1);
report = surfsmooth3d.edgepreserve.validity_report(partial,'adaptive partial surface',reference);
reportSmooth = surfsmooth3d.edgepreserve.validity_report(Ssmooth,'smooth companion (not tolerance-certified)',reference);
beta = cellfun(@(x) x.beta,records,'UniformOutput',false);
sigma = cellfun(@(x) x.sigma,records,'UniformOutput',false);
result = struct('surface',partial,'smooth',Ssmooth,'reference',reference,'report',report, ...
    'report_smooth',reportSmooth,'beta',[beta{:}].','sigma',[sigma{:}].', ...
    'info',info,'settings',s);
if ~info.converged
    warning('adaptiveblend:toleranceUnmet','%s: %d partial-surface patches unresolved.',reason,sum(unresolved));
end
end

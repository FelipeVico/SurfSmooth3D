function manifest = save_result(pkg,result,outputParent)
%SAVE_RESULT Publish each geometrically valid member; tolerance is for partial only.
if ~result.report.is_valid && ~result.report_smooth.is_valid
    error('adaptiveblend:invalidGeometry','Neither member of the pair is geometrically valid; export refused.');
end
s = result.settings; info = result.info;
parent = fullfile(outputParent,pkg.geometry_name);
if ~isfolder(parent), mkdir(parent); end
staging = tempname(parent); mkdir(staging);
cleanup = onCleanup(@() remove_staging(staging));
number = @(x) strrep(strrep(sprintf('%.8g',x),'.','p'),'-','m');
ids = 'none';
if ~isempty(s.selectedEdgeIds), ids = strjoin(arrayfun(@num2str,unique(s.selectedEdgeIds),'UniformOutput',false),'-'); end
% Bound filename length; the complete selection is always saved in parameters.mat.
if numel(ids)>60, ids = sprintf('%dedges',numel(s.selectedEdgeIds)); end
base = sprintf('%s_twostage_adaptiveblend_p%d_rlam%s_adapt%d_tol%s', ...
    pkg.name,s.order,number(s.rlam),s.adapt_sigma,number(s.eps_adapt));
suffix = ''; if ~info.converged, suffix = '_tolunmet'; end
names = {sprintf('%s_edgepreserve_edges%s_ceps%s_crad%s_ctau%s%s.go3', ...
    base,ids,number(s.c_eps),number(s.c_rad),number(s.c_tau),suffix); ...
    [base '_fullysmooth_companion.go3']};
surfaces = {result.surface,result.smooth}; reports = {result.report,result.report_smooth};
status = strings(2,1); messages = strings(2,1);
for k = 1:2
    if ~reports{k}.is_valid
        status(k) = "rejected"; messages(k) = "Geometry validity gate failed."; continue;
    end
    try
        surfsmooth3d.edgepreserve.export_surface(fullfile(staging,names{k}),surfaces{k},reports{k});
        status(k) = "saved";
    catch exception
        status(k) = "failed"; messages(k) = string(exception.message);
        candidate = fullfile(staging,names{k});
        if isfile(candidate), delete(candidate); end
    end
end
if ~any(status=="saved"), error('adaptiveblend:export','No surface could be exported: %s',strjoin(messages,'; ')); end
manifest = table(["partial";"fully smooth companion"],string(names),status,messages, ...
    'VariableNames',{'surface','filename','status','message'});
parameters = struct('source_paths',{pkg.source_paths},'source_hashes',{pkg.source_hashes}, ...
    'source_order',pkg.source_order,'source_level',pkg.source_level,'drawn_settings',s, ...
    'sigma_mode',surfsmooth3d.edgepreserve.sigma_mode_label(s.adapt_sigma),'info',info, ...
    'partial_validity',result.report.stats,'smooth_validity',result.report_smooth.stats, ...
    'tolerance_applies_to','partial only; smooth companion not tolerance-certified', ...
    'independent_refinement',true,'tolerance_is_bie_error_bound',false);
save(fullfile(staging,'parameters.mat'),'parameters');
leafTable = array2table([info.parent_ids info.depth reshape(info.maps,6,[]).' ...
    info.indicators info.launch_diameter info.unresolved info.independent_checked], ...
    'VariableNames',{'coarse_parent','depth','u0','v0','uu','vu','uv','vv', ...
    'position_tail','normal_tail','position_check','normal_check','launch_diameter','unresolved','checked'});
writetable(leafTable,fullfile(staging,'leaf_diagnostics.csv'));
writetable(info.history,fullfile(staging,'convergence_history.csv'));
writetable(manifest,fullfile(staging,'manifest.csv'));
[~,tag] = fileparts(staging);
destination = fullfile(parent,['adaptive_blend_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '_' tag]);
[ok,message] = movefile(staging,destination);
if ~ok, error('adaptiveblend:publish','%s',message); end
manifest.filename = fullfile(string(destination),manifest.filename);
fprintf('Saved adaptive blend pair: %s\n',destination);
if any(status~="saved"), warning('adaptiveblend:partialExport','Some pair members were not exported; see manifest.csv.'); end
end

function remove_staging(folder)
if isfolder(folder), rmdir(folder,'s'); end
end

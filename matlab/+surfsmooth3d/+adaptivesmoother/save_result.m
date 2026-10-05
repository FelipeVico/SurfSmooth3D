function filename = save_result(pkg, result, outputParent)
%SAVE_RESULT Validate and publish a complete adaptive run, without overwriting.
if ~result.report.is_valid
    error('adaptivesmoother:invalidSurface','Refusing export: %d invalid patches.', ...
        result.report.stats.bad_patch_count);
end
settings = result.settings;
info = result.info;
suffix = '';
if ~info.converged, suffix = '_tolunmet'; end
number = @(x) strrep(strrep(sprintf('%.8g',x),'.','p'),'-','m');
name = sprintf('%s_fullysmooth_twostage_adaptive_p%d_rlam%s_adapt%d_tol%s%s', ...
    pkg.name,settings.order,number(settings.rlam),settings.adapt_sigma, ...
    number(settings.eps_adapt),suffix);
parent = fullfile(outputParent,pkg.geometry_name);
if ~isfolder(parent), mkdir(parent); end
staging = tempname(parent);
mkdir(staging);
cleanup = onCleanup(@() remove_staging(staging));
parameters = struct('source_paths',{pkg.source_paths},'source_hashes',{pkg.source_hashes}, ...
    'source_order',pkg.source_order,'source_level',pkg.source_level, ...
    'drawn_settings',settings,'sigma_mode',surfsmooth3d.edgepreserve.sigma_mode_label(settings.adapt_sigma), ...
    'info',info,'validity',result.report.stats, ...
    'independent_refinement',true,'tolerance_is_bie_error_bound',false);
save(fullfile(staging,'parameters.mat'),'parameters');
leafTable = array2table([info.parent_ids,info.depth,reshape(info.maps,6,[]).', ...
    info.indicators,info.launch_diameter,info.unresolved,info.independent_checked], ...
    'VariableNames',{'coarse_parent','depth','map_u0','map_v0','map_uu','map_vu', ...
    'map_uv','map_vv','position_tail','normal_tail','position_check','normal_check', ...
    'launch_diameter','unresolved','independent_checked'});
writetable(leafTable,fullfile(staging,'leaf_diagnostics.csv'));
writetable(info.history,fullfile(staging,'convergence_history.csv'));
surfsmooth3d.edgepreserve.export_surface(fullfile(staging,[name '.go3']),result.surface,result.report);
[~,uniqueName] = fileparts(staging);
timestamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
destination = fullfile(parent,['adaptive_' timestamp '_' uniqueName]);
[ok,message] = movefile(staging,destination);
if ~ok, error('adaptivesmoother:publish','%s',message); end
filename = fullfile(destination,[name '.go3']);
fprintf('Saved adaptive surface: %s\n',filename);
end

function remove_staging(directory)
if isfolder(directory), rmdir(directory,'s'); end
end

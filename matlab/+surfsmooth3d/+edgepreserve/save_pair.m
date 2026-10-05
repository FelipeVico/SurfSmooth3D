function rows = save_pair(runDir, pkg, result)
%SAVE_PAIR Validate and save full and partial surfaces independently.
rows = repmat(empty_row(),2,1);
modes = {'fullysmooth','edgepreserve'};
surfaces = {result.S_smooth,result.S_blend};
reports = {result.report_smooth,result.report_blend};
for k = 1:2
    row = empty_row();
    row.order = result.order;
    row.refinement = result.level;
    row.mode = string(modes{k});
    row.filename = string(surfsmooth3d.edgepreserve.surface_filename( ...
        pkg.geometry_name,result.params,result.order,result.level,modes{k}));
    report = reports{k};
    if ~isempty(report)
        row.normal_cross_dot_min = report.stats.normal_cross_dot_min;
        row.jac_ratio_min = report.stats.jac_ratio_min;
        row.max_abs_curvature = report.stats.mean_curv_abs_max;
        row.bad_nodes = report.stats.bad_node_count;
        row.bad_patches = report.stats.bad_patch_count;
    end
    try
        if isempty(surfaces{k})
            error('edgepreserve:blendBuild', '%s',result.blend_error);
        end
        if ~report.is_valid || report.has_warnings
            bad = surfsmooth3d.edgepreserve.bad_nodes_table(report);
            writetable(bad,fullfile(runDir, ...
                replace(char(row.filename),'.go3','_bad_nodes.csv')));
        end
        surfsmooth3d.edgepreserve.export_surface(fullfile(runDir,row.filename),surfaces{k},report);
        row.status = "saved";
        row.validity = "PASS";
        if report.has_warnings, row.validity = "WARN"; end
    catch exception
        row.status = "failed";
        row.validity = "FAIL";
        row.error_message = string(exception.message);
    end
    rows(k) = row;
end
end

function row = empty_row()
row = struct('order',0,'refinement',0,'mode',"",'filename',"", ...
    'status',"failed",'validity',"UNAVAILABLE",'normal_cross_dot_min',nan, ...
    'jac_ratio_min',nan,'max_abs_curvature',nan,'bad_nodes',nan, ...
    'bad_patches',nan,'error_message',"");
end

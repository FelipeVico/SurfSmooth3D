function result = update_blend(result, edges, settings)
%UPDATE_BLEND Change cheap controls without changing the drawn smoother state.
for name = {'selectedEdgeIds','c_eps','c_rad','c_tau'}
    result.params.(name{1}) = settings.(name{1});
end
result.S_blend = [];
result.report_blend = [];
result.beta = [];
result.dsoft = [];
result.blend_error = '';
try
    if isempty(result.sigma) && ~isempty(settings.selectedEdgeIds)
        error('edgepreserve:missingSigma', ...
            'Sigma evaluation failed. Draw / Recompute before blending selected edges.');
    end
    [result.S_blend,result.beta,result.dsoft] = surfsmooth3d.edgepreserve.blend( ...
        result.S_cad,result.S_smooth,edges,result.params);
    result.report_blend = surfsmooth3d.edgepreserve.validity_report( ...
        result.S_blend,'edgepreserve',result.S_cad);
catch exception
    result.blend_error = exception.message;
end
end

function filename = surface_filename(geometryName, params, order, level, mode)
%SURFACE_FILENAME Include output resolution and drawn, not pending, parameters.
geometryName = regexprep(geometryName,'[^A-Za-z0-9_-]','_');
base = sprintf('%s_p%d_rl%d_%s_twostage_rlam%s_adapt%d', ...
    geometryName,order,level,mode,number_text(params.rlam),params.adapt_sigma);
if strcmp(mode,'edgepreserve')
    ids = strjoin(arrayfun(@num2str,params.selectedEdgeIds, ...
        'UniformOutput',false),'-');
    if isempty(ids), ids = 'none'; end
    base = sprintf('%s_edges%s_ceps%s_crad%s_ctau%s',base,ids, ...
        number_text(params.c_eps),number_text(params.c_rad),number_text(params.c_tau));
end
filename = [base '.go3'];
end

function text = number_text(value)
text = sprintf('%.6g',value);
text = strrep(strrep(strrep(text,'-','m'),'+',''),'.','p');
end

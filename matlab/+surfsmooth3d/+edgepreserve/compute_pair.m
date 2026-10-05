function result = compute_pair(pkg, work, S_smooth, level, settings, sigmaFcn)
%COMPUTE_PAIR Align the fixed CAD source, evaluate sigma, and form both outputs.
if nargin < 6, sigmaFcn = @surfsmooth3d.multiscale_mesher_sigma_eval; end
S_cad = surfsmooth3d.edgepreserve.resample_cad(pkg,S_smooth,level);
opts = surfsmooth3d.edgepreserve.mesher_options(pkg,work,settings,level);
params = settings;
params.sigma = [];
params.geometryDiameter = max(norm(max(pkg.surface.r,[],2)- ...
    min(pkg.surface.r,[],2)),1);
result = struct('S_cad',S_cad,'S_smooth',S_smooth,'S_blend',[], ...
    'beta',[],'dsoft',[],'sigma',[],'gradSigma',[], ...
    'params',params,'order',double(S_smooth.norders(1)),'level',level, ...
    'report_smooth',surfsmooth3d.edgepreserve.validity_report(S_smooth,'fullysmooth',S_cad), ...
    'report_blend',[],'blend_error','');
try
    [sigma,gradSigma] = sigmaFcn(work.scaffold_file,pkg.source_order,opts,S_cad.r);
    if numel(sigma) ~= S_cad.npts || any(~isfinite(sigma(:)) | sigma(:)<=0) || ...
            ~isequal(size(gradSigma),size(S_cad.r)) || any(~isfinite(gradSigma(:)))
        error('edgepreserve:invalidSigma', 'Sigma evaluator returned invalid data.');
    end
    result.sigma = sigma(:);
    result.gradSigma = gradSigma;
    result.params.sigma = sigma(:);
    result = surfsmooth3d.edgepreserve.update_blend(result,pkg.edges,settings);
catch exception
    % A failed partial-field evaluation must not discard a valid full surface.
    result.blend_error = exception.message;
    if isempty(settings.selectedEdgeIds)
        result = surfsmooth3d.edgepreserve.update_blend(result,pkg.edges,settings);
    end
end
end

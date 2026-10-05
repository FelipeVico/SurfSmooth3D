function [sigma, gradSigma, meta] = multiscale_mesher_sigma_eval(fnamein, norder, opts, targets)
%MULTISCALE_MESHER_SIGMA_EVAL Evaluate surface-smoother sigma at targets.
%
% Syntax
%   [sigma, gradSigma] = surfsmooth3d.multiscale_mesher_sigma_eval(fnamein, norder, opts, targets)
%
% targets must be a 3 x ntargets array. The geometry and options are
% interpreted in the same way as multiscale_mesher.
% opts.adapt_sigma: 0 constant mean longest-side sigma, 1 longest-side
% adaptive (default), 2 recursive adaptive, or 3 shortest-side adaptive.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if nargin < 4
        error('MULTISCALE_MESHER_SIGMA_EVAL: missing target array');
    end
    if size(targets, 1) ~= 3
        error('MULTISCALE_MESHER_SIGMA_EVAL: targets must be 3 x ntargets');
    end

    d = dir(fnamein);
    if isempty(d)
        error('MULTISCALE_MESHER_SIGMA_EVAL: invalid file');
    end
    fnameuse = fullfile(d.folder, d.name);
    targets = double(targets);
    ntarg = size(targets, 2);

    norder_skel = 12;
    if isfield(opts, 'nquad')
        norder_skel = opts.nquad;
    end
    norder_smooth = norder;

    adapt_sigma = 1;
    if isfield(opts, 'adapt_sigma')
        adapt_sigma = opts.adapt_sigma;
    end
    if (~isnumeric(adapt_sigma) && ~islogical(adapt_sigma)) || ...
            ~isscalar(adapt_sigma) || ~isreal(adapt_sigma) || ~ismember(adapt_sigma, 0:3)
        error('MULTISCALE_MESHER_SIGMA_EVAL:InvalidSigmaMode', ...
            'opts.adapt_sigma must be 0, 1, 2, or 3.');
    end
    adapt_sigma = double(adapt_sigma);

    rlam = 10;
    if isfield(opts, 'rlam')
        rlam = opts.rlam;
    end

    if isfield(opts, 'filetype')
        ifiletype = opts.filetype;
    else
        ier = 0;
        ifiletype = 0;
        mex_id_ = 'MWF77_get_filetype(c i cstring[x], c io int64_t[x], c io int64_t[x])';
[ifiletype, ier] = surfsmooth3d_routs(mex_id_, fnameuse, ifiletype, ier, 1000, 1, 1);
        if ier > 0
            error('MULTISCALE_MESHER_SIGMA_EVAL: error determining file type');
        end
    end

    if norder > 20 || norder < 1
        error('MULTISCALE_MESHER_SIGMA_EVAL: norder must be between 1 and 20');
    end
    if norder_skel > 20 || norder_skel < 1
        error('MULTISCALE_MESHER_SIGMA_EVAL: opts.nquad must be between 1 and 20');
    end

    ifcad = 0;
    fcad = 'tmp';
    if isfield(opts, 'fcad')
        fcad = opts.fcad;
        ifcad = 1;
    end

    sigma = zeros(ntarg, 1);
    gradSigma = zeros(3, ntarg);
    ier = 0;

    mex_id_ = 'MWF77_multiscale_mesher_sigma_eval(c i cstring[x], c i int64_t[x], c i int64_t[x], c i cstring[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i double[xx], c io double[x], c io double[xx], c io int64_t[x])';
[sigma, gradSigma, ier] = surfsmooth3d_routs(mex_id_, fnameuse, ifiletype, ifcad, fcad, norder_skel, norder_smooth, adapt_sigma, rlam, ntarg, targets, sigma, gradSigma, ier, 1000, 1, 1, 1000, 1, 1, 1, 1, 1, 3, ntarg, ntarg, 3, ntarg, 1);

    if ier > 0
        error('MULTISCALE_MESHER_SIGMA_EVAL: Fortran sigma evaluation failed with ier=%d', ier);
    end

    meta = struct();
    meta.filetype = ifiletype;
    meta.ifcad = ifcad;
    meta.norder_skel = norder_skel;
    meta.norder_smooth = norder_smooth;
    meta.adapt_sigma = adapt_sigma;
    meta.rlam = rlam;
end


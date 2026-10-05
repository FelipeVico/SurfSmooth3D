function S = multiscale_mesher(fnamein, norder, opts)
% 
%  MULTISCALE_MESHER creates a smooth high order triangulated surface 
%    based on an input mesh file containing first/second order
%    triangles/quads. If the mesh comprises of quads, then each quad
%    is split into two triangles. 
%
%  Supported mesh formats include, .gidmsh, .tri, .msh, gmshv2, 
%  and gmshv4
%  
%  Syntax
%    S = surfsmooth3d.multiscale_mesher(fnamein, norder)
%    S = surfsmooth3d.multiscale_mesher(fnamein, norder, opts)
%
%  Input arguments:
%    * fnamein: input mesh file name
%    * norder: order of discretization for smooth surface
%    * opts: options struct (optional)
%        opts.nquad (12), quadrature order for computing level set
%        opts.rlam (10), smoothing parameter (should be between 2.5 and 10)
%        opts.adapt_sigma (1), 0: constant mean longest-side sigma
%                              1: longest-side adaptive sigma (default)
%                              2: existing recursive adaptive sigma
%                              3: shortest-side adaptive sigma
%        opts.nrefine (0), number of refinements
%        opts.fcad (), cad file info to be used for defining surface for 
%        computing level set
%        opts.two_stage_smoother (false), first project scaffold vertices
%        using the CAD-skeleton level set, then run the full smoother
%        opts.filetype, type of file
%          filetype = 1, for .msh from gidmsh
%          filetype = 2, for .tri
%          filetype = 3, for .gidmsh
%          filetype = 4, for .msh gmsh v2
%          filetype = 5, for .msh gmsh v2
%                               
%      
%  
    d = dir(fnamein);
    if isempty(d)
        error('MULTISCALE_MESHER: invalid file\n');
    end
    fnameuse = fullfile(d.folder, d.name);

    if nargin < 3
        opts = [];
    end


    fnameoutuse = fullfile(d.folder, '/tmp');
    
    norder_skel = 12;
    if isfield(opts, 'nquad')
        norder_skel = opts.nquad;
    end

    norder_smooth = norder;
    nrefine = 0;
    if isfield(opts, 'nrefine')
        nrefine = opts.nrefine;
    end


    adapt_sigma = 1;
    if isfield(opts, 'adapt_sigma')
        adapt_sigma = opts.adapt_sigma;
    end
    if (~isnumeric(adapt_sigma) && ~islogical(adapt_sigma)) || ...
            ~isscalar(adapt_sigma) || ~isreal(adapt_sigma) || ~ismember(adapt_sigma, 0:3)
        error('MULTISCALE_MESHER:InvalidSigmaMode', ...
            'opts.adapt_sigma must be 0, 1, 2, or 3.');
    end
    adapt_sigma = double(adapt_sigma);
    rlam = 10;
    if isfield(opts, 'rlam')
        rlam = opts.rlam;
    end
    two_stage_smoother = false;
    if isfield(opts, 'two_stage_smoother')
        two_stage_smoother = logical(opts.two_stage_smoother);
    end

    if isfield(opts, 'filetype')
        ifiletype = opts.filetype;
    else
        ier = 0;
        ifiletype = 0;
        mex_id_ = 'MWF77_get_filetype(c i cstring[x], c io int64_t[x], c io int64_t[x])';
[ifiletype, ier] = surfsmooth3d_routs(mex_id_, fnameuse, ifiletype, ier, 1000, 1, 1);
        if ier > 0
            error('MULTISCALE_MESHER: error determining file type\n');
        end
    end

    if norder > 20 || norder < 1
        error('MULTISCALE_MESHER: norder too high, must be less than 20');
    end
    
    if norder_skel > 20 || norder_skel < 1
        error('MULTISCALE_MESHER: opts.nquad too high, must be less than 20');
    end
    ier = 0;
    ifcad = 0;
    fcad = 'tmp';
    if isfield(opts, 'fcad')
      fcad = opts.fcad;
      ifcad = 1;
    end
    if two_stage_smoother && ifcad == 0
        error(['MULTISCALE_MESHER: opts.two_stage_smoother requires ', ...
            'opts.fcad for the CAD skeleton level-set source']);
    end
    reportFile = [fnameoutuse '_newton_failure.txt'];
    if two_stage_smoother
        if isfile(reportFile), delete(reportFile); end
        if isfile(reportFile)
            error('MULTISCALE_MESHER:NewtonDiagnosticFile', ...
                'Cannot remove a stale Newton diagnostic: %s',reportFile);
        end
        mex_id_ = 'MWF77_multiscale_mesher_two_stage(c i cstring[x], c i int64_t[x], c i int64_t[x], c i cstring[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i cstring[x], c io int64_t[x])';
[ier] = surfsmooth3d_routs(mex_id_, fnameuse, ifiletype, ifcad, fcad, norder_skel, norder_smooth, nrefine, adapt_sigma, rlam, fnameoutuse, ier, 1000, 1, 1, 1000, 1, 1, 1, 1, 1, 1000, 1);
    else
        mex_id_ = 'MWF77_multiscale_mesher(c i cstring[x], c i int64_t[x], c i int64_t[x], c i cstring[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i cstring[x], c io int64_t[x])';
[ier] = surfsmooth3d_routs(mex_id_, fnameuse, ifiletype, ifcad, fcad, norder_skel, norder_smooth, nrefine, adapt_sigma, rlam, fnameoutuse, ier, 1000, 1, 1, 1000, 1, 1, 1, 1, 1, 1000, 1);
    end

    if ier == 7
        exception = MException('MULTISCALE_MESHER:NewtonRadiusGuard', ...
            ['Two-stage Newton stopped: a target left the fixed 10R sphere ', ...
             'or became nonfinite. No failed surface was returned. Diagnostic: %s'],reportFile);
        exception = addCause(exception,MException( ...
            'MULTISCALE_MESHER:NewtonDiagnosticFile','%s',reportFile));
        throw(exception);
    end
    if ier == 3 || ier == 4
        if two_stage_smoother
            reason = newton_failure_reason(reportFile,ier);
            exception = MException('MULTISCALE_MESHER:NewtonProjectionFailure', ...
                ['Two-stage Newton stopped: %s. No failed surface was returned.\n', ...
                 'A nearby Phi=1/2 intersection along a pseudonormal is not guaranteed. ', ...
                 'Increasing rlam reduces sigma; recompute and check geometry validity.\n', ...
                 'Diagnostic: %s'],reason,reportFile);
            exception = addCause(exception,MException( ...
                'MULTISCALE_MESHER:NewtonDiagnosticFile','%s',reportFile));
            throw(exception);
        end
        error_message = ['MULTISCALE_MESHER: error in main smoothing routine', ...
                'try a larger value of rlam (if not set, use opts.rlam and set it to ',...
                'a value greater than 10)\n' ...
                'if that does not work, then mesh cannot be smooth with the surface smoother\n'];
        error(error_message);
    end
    if ier == 6
        error(['MULTISCALE_MESHER: two-stage vertex projection produced ', ...
            'a degenerate scaffold']);
    end

    S = cell(1,nrefine+1); 
    for i=0:nrefine
        fname = fullfile(d.folder, ['/tmp_o' num2str(norder,'%02.f') '_r' num2str(i,'%02.f') '.go3']);
        S{i+1} = surfsmooth3d.surfer.load_from_file(fname);
        delete(fname);
    end
end

function reason = newton_failure_reason(filename,ier)
reason = 'Newton failed to converge';
if ier == 4, reason = 'Newton correction or directional derivative failed its stopping criterion'; end
fid = fopen(filename,'r');
if fid < 0, return; end
cleanup = onCleanup(@() fclose(fid));
if strcmp(fgetl(fid),'NEWTON_PROJECTION_FAILURE_V2')
    line = fgetl(fid);
    prefix = 'failure_reason ';
    if ischar(line) && startsWith(line,prefix)
        reason = strtrim(line(numel(prefix)+1:end));
    end
end
end

function [S, info] = multiscale_mesher_adaptive(scaffoldFile, norder, opts)
%MULTISCALE_MESHER_ADAPTIVE Fixed-order adaptive two-stage CAD smoothing.
% [S,info] = surfsmooth3d.multiscale_mesher_adaptive(scaffoldFile,norder,opts)
% opts.fcad is required (x y z nx ny nz w point-skeleton file).
% Existing options: nquad=12, rlam=10, adapt_sigma=1, filetype inferred.
% Sigma modes: 0 constant, 1 longest-side, 2 recursive, 3 shortest-side.
% Adaptive options: eps_adapt=1e-4, max_refine=3, max_points=2000000.
% Order is 1:20, maximum depth is 0:20. Leaves refine independently.
% The CAD source and original scaffold sigma field stay fixed. Tolerance
% measures geometric resolution, not BIE error. Nonconvergence at a resource
% limit returns the last completed surface with info.converged=false.
% info.maps(:,1,k)+info.maps(:,2:3,k)*uv maps leaf k to its coarse parent.
% info.indicators(:,1:4): position RMS tail/h, unit-normal RMS tail,
% independent Phi distance/h, and independent normal discrepancy.
% A negative independent indicator means not evaluated (a tail already failed).
% info.independent_checked disambiguates those entries; accepted leaves have
% all four checks. info.history records per-pass maxima, counts and depths.
% Hanging reference edges are allowed; exact interpatch continuity is not
% imposed. A geometrically valid mesh can still have info.converged=false.
if nargin < 3, opts = struct(); end
if ~isstruct(opts) || ~isscalar(opts) || ~isfield(opts,'fcad') || ~isfile(opts.fcad)
    error('adaptivesmoother:cadSource','An existing opts.fcad point skeleton is required.');
end
if ~isfile(scaffoldFile)
    error('adaptivesmoother:scaffold','Scaffold file does not exist.');
end
defaults = struct('nquad',12,'rlam',10,'adapt_sigma',1, ...
    'eps_adapt',1e-4,'max_refine',3,'max_points',2000000);
fields = fieldnames(defaults);
for k = 1:numel(fields)
    if ~isfield(opts,fields{k}), opts.(fields{k}) = defaults.(fields{k}); end
end
validateattributes(norder,{'numeric'},{'real','scalar','integer','>=',1,'<=',20});
validateattributes(opts.nquad,{'numeric'},{'real','scalar','integer','>=',1,'<=',20});
validateattributes(opts.rlam,{'numeric'},{'real','scalar','finite','positive'});
validateattributes(opts.adapt_sigma,{'numeric','logical'},{'real','scalar','integer','>=',0,'<=',3});
validateattributes(opts.eps_adapt,{'numeric'},{'real','scalar','finite','>',0,'<',1});
validateattributes(opts.max_refine,{'numeric'},{'real','scalar','integer','>=',0,'<=',20});
validateattributes(opts.max_points,{'numeric'},{'real','scalar','integer','>=',1,'<=',flintmax});
if ~isfield(opts,'filetype')
    id = 'MWF77_get_filetype(c i cstring[x], c io int64_t[x], c io int64_t[x])';
    [opts.filetype,ier] = surfsmooth3d_routs(id,char(scaffoldFile),0,0,1000,1,1);
    if ier ~= 0, error('adaptivesmoother:filetype','Cannot determine scaffold type.'); end
end
validateattributes(opts.filetype,{'numeric'},{'real','scalar','integer','>=',1,'<=',5});
if exist('surfsmooth3d_adaptive_routs','file') ~= 3
    error('adaptivesmoother:build','Build the new gateway with make matlab first.');
end
directory = tempname;
mkdir(directory);
cleanup = onCleanup(@() rmdir(directory,'s'));
root = fullfile(directory,'adaptive');
ier = surfsmooth3d_adaptive_routs(char(scaffoldFile),char(opts.fcad),root, ...
    double(opts.filetype),double(opts.nquad),double(norder),double(opts.adapt_sigma), ...
    double(opts.rlam),double(opts.eps_adapt),double(opts.max_refine),double(opts.max_points));
if ier ~= 0
    message = sprintf('Adaptive smoother failed (code %d). No failed mesh was returned.',ier);
    identifier = 'adaptivesmoother:solve';
    if ier == 7, identifier = 'MULTISCALE_MESHER:NewtonRadiusGuard'; end
    if ismember(ier,[3 4]), identifier = 'MULTISCALE_MESHER:NewtonProjectionFailure'; end
    if ier == 6, message = 'Vertex projection produced a degenerate scaffold.'; end
    if ier == 8, message = 'Level-zero mesh already exceeds opts.max_points.'; end
    if ier == 9, message = 'Cannot write adaptive result files.'; end
    exception = MException(identifier,'%s',message);
    reportFile = [root '_newton_failure.txt'];
    if isfile(reportFile)
        preserved = [tempname '_newton_failure.txt'];
        copyfile(reportFile,preserved);
        exception = addCause(exception,MException( ...
            'MULTISCALE_MESHER:NewtonDiagnosticFile','%s',preserved));
    end
    throw(exception);
end
S = surfsmooth3d.surfer.load_from_file([root '.go3']);
leaves = load([root '_leaves.txt'],'-ascii');
history = load([root '_history.txt'],'-ascii');
status = load([root '_status.txt'],'-ascii');
if ~ismember(status(1),0:2) || size(leaves,1) ~= S.npatches || size(leaves,2) ~= 13
    error('adaptivesmoother:output','Malformed adaptive diagnostics.');
end
reasons = {'tolerance achieved','maximum depth reached','point budget reached'};
info = struct('converged',status(1)==0,'stop_reason',reasons{status(1)+1}, ...
    'parent_ids',leaves(:,1),'depth',leaves(:,2), ...
    'maps',reshape(leaves(:,3:8).',2,3,[]),'indicators',leaves(:,9:12), ...
    'launch_diameter',leaves(:,13),'center',status(2:4).','radius',status(5), ...
    'order',norder,'settings',opts);
info.indicator_names = {'position_tail','normal_tail','position_check','normal_check'};
info.independent_checked = all(info.indicators(:,3:4)>=0,2);
info.unresolved = max(info.indicators,[],2)>opts.eps_adapt;
info.unresolved_count = nnz(info.unresolved);
info.achieved_depth = max(info.depth);
info.history = array2table(history,'VariableNames',{'pass','patches','points', ...
    'marked','unresolved','depth','position_tail','normal_tail','position_check', ...
    'normal_check','independent_checked'});
if ~info.converged
    warning('adaptivesmoother:toleranceUnmet','%s: %d unresolved patches.', ...
        info.stop_reason,info.unresolved_count);
end
end

function [manifest, runDir] = run_batch(pkg, settings, outputParent, progress, hooks)
%RUN_BATCH Solve once per order, keeping CAD/scaffold/edges fixed throughout.
% Optional function handles support small deterministic failure-path tests.
if nargin < 4 || isempty(progress), progress = @(text) fprintf('%s\n',text); end
if nargin < 5, hooks = struct(); end
if ~isfield(hooks,'solve'), hooks.solve = @surfsmooth3d.multiscale_mesher; end
if ~isfield(hooks,'sigma'), hooks.sigma = @surfsmooth3d.multiscale_mesher_sigma_eval; end
validateattributes(settings.orders,{'numeric'},{'vector','nonempty','integer','>=',1,'<=',20});
validateattributes(settings.levels,{'numeric'},{'vector','nonempty','integer','nonnegative'});
settings.orders = unique(settings.orders);
settings.levels = unique(settings.levels);
runDir = surfsmooth3d.edgepreserve.new_run_directory(outputParent,pkg,settings,'batch');
[work,cleanup] = surfsmooth3d.edgepreserve.stage_source(pkg); %#ok<ASGLU>
rows = [];
maxLevel = max(settings.levels);
for order = settings.orders
    progress(sprintf('Solving order %d through refinement %d...',order,maxLevel));
    opts = surfsmooth3d.edgepreserve.mesher_options(pkg,work,settings,maxLevel);
    solveError = '';
    try
        surfaces = hooks.solve(work.scaffold_file,order,opts);
        if numel(surfaces) ~= maxLevel+1
            error('edgepreserve:solverLevels', 'Smoother did not return all refinement levels.');
        end
    catch exception
        preserve_failure(exception,order,maxLevel,'maximum');
        solveError = exception.message;
        surfaces = {};
        progress(sprintf('Order %d maximum-level solve failed; retrying lower requested levels.',order));
    end
    if isempty(solveError)
        unneeded = setdiff(0:maxLevel,settings.levels);
        surfaces(unneeded+1) = {[]};
    end
    for level = settings.levels
        progress(sprintf('Processing order %d, refinement %d...',order,level));
        try
            if ~isempty(solveError)
                if level == maxLevel
                    error('edgepreserve:maximumSolve', '%s',solveError);
                end
                opts.nrefine = level;
                retry = hooks.solve(work.scaffold_file,order,opts);
                S = retry{level+1};
                clear retry;
            else
                S = surfaces{level+1};
                surfaces{level+1} = [];
            end
            result = surfsmooth3d.edgepreserve.compute_pair(pkg,work,S,level,settings,hooks.sigma);
            pair = surfsmooth3d.edgepreserve.save_pair(runDir,pkg,result);
            clear S result;
        catch exception
            preserve_failure(exception,order,level,'retry');
            pair = failed_pair(pkg,settings,order,level,exception.message);
        end
        rows = [rows;pair]; %#ok<AGROW>
        manifest = struct2table(rows);
        writetable(manifest,fullfile(runDir,'batch_manifest.csv'));
        progress(sprintf('Order %d, refinement %d: full %s; partial %s.', ...
            order,level,pair(1).status,pair(2).status));
    end
    clear surfaces;
end
manifest = struct2table(rows);
progress(sprintf('Saved %d/%d surfaces. Output: %s', ...
    nnz(manifest.status == "saved"),height(manifest),runDir));

    function preserve_failure(exception,order,level,attempt)
        filename = fullfile(runDir,sprintf('newton_failure_p%02d_r%02d_%s.txt',order,level,attempt));
        [report,filename] = surfsmooth3d.edgepreserve.capture_newton_failure(exception,filename);
        if isempty(report), return; end
        if ~isempty(filename)
            progress(sprintf('Newton failure diagnostic preserved: %s',filename));
        else
            progress('Newton diagnostic could not be archived; the current snapshot is still available.');
        end
        if isfield(hooks,'failure')
            try
                hooks.failure(report,filename);
            catch diagnosticError
                warning('edgepreserve:newtonDiagnostic','Failure-view callback failed: %s',diagnosticError.message);
            end
        end
    end
end

function rows = failed_pair(pkg,settings,order,level,message)
modes = {'fullysmooth','edgepreserve'};
for k = 1:2
    row = struct('order',order,'refinement',level,'mode',string(modes{k}), ...
        'filename',string(surfsmooth3d.edgepreserve.surface_filename( ...
            pkg.geometry_name,settings,order,level,modes{k})), ...
        'status',"failed",'validity',"UNAVAILABLE",'normal_cross_dot_min',nan, ...
        'jac_ratio_min',nan,'max_abs_curvature',nan,'bad_nodes',nan, ...
        'bad_patches',nan,'error_message',string(message));
    rows(k,1) = row; %#ok<AGROW>
end
end

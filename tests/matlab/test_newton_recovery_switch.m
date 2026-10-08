function test_newton_recovery_switch()
%TEST_NEWTON_RECOVERY_SWITCH Real disappearing-levelset failures through all APIs.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
package = fullfile(root,'tests','fixtures','cad', ...
    'convergence_intersection_two_balls','p4_mf0p10000_rl0');
pkg = surfsmooth3d.edgepreserve.load_package(package);
[work,cleanup] = surfsmooth3d.edgepreserve.stage_source(pkg); %#ok<ASGLU>
opts = struct('fcad',work.cad_file,'filetype',3,'nquad',pkg.source_order, ...
    'rlam',.05,'adapt_sigma',1,'two_stage_smoother',true,'nrefine',0, ...
    'max_refine',0,'eps_adapt',.01,'max_points',2000000);
lockBaseline = adaptive_blend_lock_baseline();
settings = surfsmooth3d.adaptiveblend.default_settings();
settings.order = 4; settings.rlam = .05; settings.max_refine = 0;
for enabled = [false true]
    opts.newton_recovery = enabled;
    settings.newton_recovery = enabled;
    calls = {@() surfsmooth3d.multiscale_mesher(work.scaffold_file,4,opts), ...
        @() surfsmooth3d.multiscale_mesher_adaptive(work.scaffold_file,4,opts), ...
        @() surfsmooth3d.adaptiveblend.solve(pkg,settings)};
    for k = 1:numel(calls)
        exception = expected_radius_failure(calls{k});
        report = cause_path(exception,'MULTISCALE_MESHER:NewtonRecoveryFile');
        if enabled
            assert(~isempty(report) && isfile(report),'Recovery report was not preserved.');
            info = surfsmooth3d.read_newton_recovery(report,true);
            assert(info.deferred>0 && info.recovered==0 && info.unresolved==info.deferred);
            assert(all(isfinite(vertcat(info.events(1).targets.base_point)),'all'));
            assert(all(strcmp({info.events(1).targets.outcome},'no_downward_bracket')));
            delete(report);
        else
            assert(isempty(report),'Disabled recovery still produced a recovery report.');
            if k==1, assert(~isfile(fullfile(work.directory,'tmp_newton_recovery.txt'))); end
        end
        snapshot = cause_path(exception,'MULTISCALE_MESHER:NewtonDiagnosticFile');
        assert(~isempty(snapshot) && isfile(snapshot));
        delete(snapshot);
    end
end
assert(isempty(dir(fullfile(work.directory,'*.go3'))), ...
    'A failed uniform solve exported a surface.');
assert(mislocked('surfsmooth3d_adaptive_blend_routs')==lockBaseline, ...
    'Failed blend open left a session locked.');
fprintf('PASS: enabled/disabled real 10R failures through uniform/adaptive/blend APIs; reports preserved and no failed export.\n');
end

function exception = expected_radius_failure(call)
try
    call();
catch exception
    assert(strcmp(exception.identifier,'MULTISCALE_MESHER:NewtonRadiusGuard'), ...
        'Wrong failure identifier: %s',exception.identifier);
    return
end
error('recovery_test:expectedFailure','Expected the oversmoothed level set to disappear.');
end

function filename = cause_path(exception,identifier)
filename = '';
for k = 1:numel(exception.cause)
    if strcmp(exception.cause{k}.identifier,identifier)
        filename = exception.cause{k}.message;
        return
    end
    filename = cause_path(exception.cause{k},identifier);
    if ~isempty(filename), return; end
end
end

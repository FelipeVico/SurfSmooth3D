function test_newton_guard_batch(packageDir, snapshotFile, outputParent)
%TEST_NEWTON_GUARD_BATCH Preserve each controlled failure before batch retries.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
pkg = surfsmooth3d.edgepreserve.load_package(packageDir);
settings = surfsmooth3d.edgepreserve.default_settings();
settings.orders = [4 6];
settings.levels = [0 1 2];
callbacks = {};
sourceReport = surfsmooth3d.edgepreserve.read_newton_failure(snapshotFile);
identifier = 'MULTISCALE_MESHER:NewtonRadiusGuard';
if strcmp(sourceReport.failure_kind,'newton'), identifier = 'MULTISCALE_MESHER:NewtonProjectionFailure'; end
expectedCount = nnz(sourceReport.rejected);
exception = MException(identifier,'Synthetic Newton failure.');
exception = addCause(exception,MException('MULTISCALE_MESHER:NewtonDiagnosticFile','%s',snapshotFile));
hooks = struct('solve',@fail_solve,'failure',@capture_callback);
[manifest,directory] = surfsmooth3d.edgepreserve.run_batch(pkg,settings,outputParent,[],hooks);
assert(height(manifest)==12 && all(manifest.status=="failed"));
assert(numel(callbacks)==6 && numel(dir(fullfile(directory,'newton_failure_*.txt')))==6);
for k = 1:numel(callbacks)
    assert(isfile(callbacks{k}));
    report = surfsmooth3d.edgepreserve.read_newton_failure(callbacks{k});
    assert(numel(report.index)==6 && nnz(report.rejected)==expectedCount);
end
assert(isempty(dir(fullfile(directory,'*.go3'))));
for k = 1:numel(pkg.source_paths)
    assert(strcmp(surfsmooth3d.edgepreserve.file_sha256(pkg.source_paths{k}),pkg.source_hashes{k}));
end
fprintf('PASS: six Newton snapshots preserved, twelve failures recorded, batch continues, no invalid export.\n');

    function result = fail_solve(varargin) %#ok<STOUT>
        throw(exception);
    end
    function capture_callback(report,filename)
        assert(nnz(report.rejected)==expectedCount && isfile(filename));
        callbacks{end+1} = filename;
    end
end

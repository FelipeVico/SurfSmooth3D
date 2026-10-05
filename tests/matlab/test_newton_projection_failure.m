function test_newton_projection_failure(reportRoot, outputDir)
%TEST_NEWTON_PROJECTION_FAILURE Diagnose Newton stops inside the radius guard.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
if nargin < 1, reportRoot = fullfile(tempdir,'test_projection_guard'); end
if nargin < 2, outputDir = tempname; end
if ~isfolder(outputDir), mkdir(outputDir); end
snapshotFile = [reportRoot '_newton_newton_failure.txt'];
report = surfsmooth3d.edgepreserve.read_newton_failure(snapshotFile);
assert(strcmp(report.failure_kind,'newton') && contains(report.failure_reason,'1.1'));
assert(isequal(report.index,(1:6).') && nnz(report.rejected)==2);
assert(all(report.reason(1:2)==3) && all(report.distance_over_radius<10));
assert(report.step_length(2)==8 && report.step_iteration(4)==2 && report.converged(4));
assert(~report.step_available(5) && isinf(report.step_length(6)));
exception = MException('MULTISCALE_MESHER:NewtonProjectionFailure','Synthetic Newton stop.');
exception = addCause(exception,MException('MULTISCALE_MESHER:NewtonDiagnosticFile','%s',snapshotFile));
destination = fullfile(outputDir,'preserved_newton_failure.txt');
[captured,filename] = surfsmooth3d.edgepreserve.capture_newton_failure(exception,destination);
assert(strcmp(filename,destination) && isequaln(captured.step_length,report.step_length));
S = fixture_sphere(4,1,[30;-20;10],4,1);
pkg = struct('surface',S,'scaffold',struct('nodes',S.r),'geometry_name','Newton fixture');
fig = surfsmooth3d.edgepreserve.plot_newton_failure(pkg,report,[],'off');
cleanup = onCleanup(@() close(fig));
red = findobj(fig,'Tag','Failed: Newton stopping criterion');
assert(isscalar(red) && isequal([red.XData;red.YData;red.ZData],report.base(:,1:2)));
assert(isempty(findobj(fig,'Tag','Rejected: outside 10R / nonfinite')));
colored = findobj(fig,'Tag','Active: latest Newton step');
assert(isscalar(colored) && isscalar(colored.XData));
assert(abs(colored.CData-log10(.002))<1e-12);
exportapp(fig,fullfile(outputDir,'newton_projection_failure.png'));
fprintf('PASS: ordinary Newton stops, full-step snapshot, preservation and distinct failure coloring.\n');
end

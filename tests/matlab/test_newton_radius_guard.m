function test_newton_radius_guard(reportRoot, outputDir)
%TEST_NEWTON_RADIUS_GUARD Read Fortran snapshots and exercise failure coloring.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
clear dir;
if nargin < 1, reportRoot = fullfile(tempdir,'test_projection_guard'); end
if nargin < 2, outputDir = tempname; end
if ~isfolder(outputDir), mkdir(outputDir); end
snapshotFile = [reportRoot '_snapshot_newton_failure.txt'];
report = surfsmooth3d.edgepreserve.read_newton_failure(snapshotFile);
assert(strcmp(report.stage,'rv_nodes') && report.refinement==2 && report.iteration==3);
assert(isequal(report.index,(1:6).') && nnz(report.rejected)==1);
assert(report.step_length(2)==8 && report.step_iteration(4)==2 && report.converged(4));
assert(~report.step_available(5) && isnan(report.step_length(5)));
assert(isinf(report.step_length(6)) && report.step_available(6));
nonfinite = surfsmooth3d.edgepreserve.read_newton_failure([reportRoot '_nonfinite_newton_failure.txt']);
assert(nnz(nonfinite.rejected)==3 && nnz(nonfinite.reason==2)==2);
assert(~any(nonfinite.step_available));

exception = MException('MULTISCALE_MESHER:NewtonRadiusGuard','Synthetic radius rejection.');
exception = addCause(exception,MException('MULTISCALE_MESHER:NewtonDiagnosticFile','%s',snapshotFile));
destination = fullfile(outputDir,'preserved_failure.txt');
[captured,filename] = surfsmooth3d.edgepreserve.capture_newton_failure(exception,destination);
assert(strcmp(filename,destination) && isequaln(captured.base,report.base));
assert(isequaln(surfsmooth3d.edgepreserve.read_newton_failure(destination).step_length,report.step_length));
assert(isempty(surfsmooth3d.edgepreserve.capture_newton_failure(MException('test:other','Other failure'))));

S = fixture_sphere(4,1,[30;-20;10],4,1);
pkg = struct('surface',S,'scaffold',struct('nodes',S.r),'geometry_name','Guard fixture');
fig = surfsmooth3d.edgepreserve.plot_newton_failure(pkg,report,[],'off');
cleanup = onCleanup(@() close_if_open(fig));
colored = findobj(fig,'Tag','Active: latest Newton step');
assert(isscalar(colored) && numel(colored.XData)==2);
assert(max(abs(colored.CData(:)-log10([2;.002])))<1e-12);
red = findobj(fig,'Tag','Rejected: outside 10R / nonfinite');
assert(isscalar(red) && isequal(red.CData,[1 0 0]));
assert(isequal([red.XData;red.YData;red.ZData],report.base(:,1)));
ax = ancestor(red,'axes');
assert(ax.ZLim(2)<30 && report.target(3,1)>ax.ZLim(2));
assert(isequaln(getappdata(fig,'edgepreserve_newton_failure'),report));
exportapp(fig,fullfile(outputDir,'newton_failure_snapshot.png'));
same = surfsmooth3d.edgepreserve.plot_newton_failure(pkg,nonfinite,fig,'off');
assert(same==fig && numel(findobj(fig,'Type','scatter'))==3);
assert(numel(findobj(fig,'Tag','Rejected: outside 10R / nonfinite').XData)==3);
close(fig);
fprintf('PASS: Fortran snapshot parsing, preservation, log coloring, red launch points and figure reuse.\n');
end

function close_if_open(fig)
if isgraphics(fig), close(fig); end
end

function smoke_newton_radius_guard_gui(packageDir, outputDir)
%SMOKE_NEWTON_RADIUS_GUARD_GUI Reproduce the multiscale failure through the GUI.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
clear dir;
if nargin < 1 || isempty(packageDir)
    error('surfsmooth3d:externalFixture', ...
        'Pass the optional Multiscale_step p8/r0 CAD package directory explicitly.');
end
if nargin < 2, outputDir = tempname; end
if ~isfolder(outputDir), mkdir(outputDir); end
fig = surfsmooth3d.edgepreserve.launch_gui(root,packageDir,'off',outputDir);
cleanup = onCleanup(@() close_figures(fig));
initial = getappdata(fig,'edgepreserve_source');
draw = findobj(fig,'Tag','draw');
draw.Callback(draw,[]);
report = getappdata(fig,'edgepreserve_newton_failure');
assert(~isempty(report) && strcmp(report.stage,'vertices') && report.iteration==3);
assert(nnz(report.rejected)==4 && numel(report.index)==2219);
assert(isempty(getappdata(fig,'edgepreserve_preview')));
assert(strcmp(findobj(fig,'Tag','save').Enable,'off'));
assert(isequal(getappdata(fig,'edgepreserve_source').surface.r,initial.surface.r));
failure = findall(groot,'Type','figure','Name','Two-stage Newton failure');
assert(isscalar(failure));
assert(numel(findobj(failure,'Tag','Rejected: outside 10R / nonfinite').XData)==4);
assert(~isempty(findobj(failure,'Tag','Active: latest Newton step')));
exportapp(failure,fullfile(outputDir,'multiscale_newton_failure.png'));
copyfile(report.filename,fullfile(outputDir,'multiscale_newton_failure.txt'));
fprintf('PASS: real multiscale GUI stops safely, four red points, full snapshot and no failed export.\n');
close_figures(fig);
end

function close_figures(fig)
if isgraphics(fig), close(fig); end
failure = findall(groot,'Type','figure','Name','Two-stage Newton failure');
delete(failure);
end

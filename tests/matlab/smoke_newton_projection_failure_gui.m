function smoke_newton_projection_failure_gui(packageDir, outputDir, testRecovery)
%SMOKE_NEWTON_PROJECTION_FAILURE_GUI Diagnose the shortest-side cavity stop.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
if nargin < 1 || isempty(packageDir)
    error('surfsmooth3d:externalFixture', ...
        'Pass the optional Multiscale_step p8/r0 CAD package directory explicitly.');
end
if nargin < 2, outputDir = tempname; end
if nargin < 3, testRecovery = false; end
if ~isfolder(outputDir), mkdir(outputDir); end
fig = surfsmooth3d.edgepreserve.launch_gui(root,packageDir,'off',outputDir);
cleanup = onCleanup(@() close_figures(fig));
source = getappdata(fig,'edgepreserve_source');
selector = findobj(fig,'Tag','adapt');
selector.Value = 2; selector.Callback(selector,[]);
draw = findobj(fig,'Tag','draw'); draw.Callback(draw,[]);
report = getappdata(fig,'edgepreserve_newton_failure');
assert(~isempty(report) && strcmp(report.failure_kind,'newton'));
assert(strcmp(report.stage,'vertices') && contains(report.failure_reason,'1.1'));
assert(nnz(report.rejected)>0 && numel(report.index)==2219);
assert(all(report.distance_over_radius<10) && all(report.reason(report.rejected)==3));
assert(isempty(getappdata(fig,'edgepreserve_preview')));
assert(strcmp(findobj(fig,'Tag','save').Enable,'off'));
assert(isequal(getappdata(fig,'edgepreserve_source').surface.r,source.surface.r));
failure = findall(groot,'Type','figure','Name','Two-stage Newton failure');
assert(isscalar(failure));
red = findobj(failure,'Tag','Failed: Newton stopping criterion');
assert(numel(red.XData)==nnz(report.rejected));
assert(~isempty(findobj(failure,'Tag','Active: latest Newton step')));
status = strjoin(cellstr(findobj(fig,'Tag','status').String),newline);
assert(contains(status,'Newton correction') && ~contains(status,'\n'));
exportapp(failure,fullfile(outputDir,'multiscale_shortest_newton_failure.png'));
copyfile(report.filename,fullfile(outputDir,'multiscale_shortest_newton_failure.txt'));
if testRecovery
    rlam = findobj(fig,'Tag','rlam');
    rlam.String = '6'; rlam.Callback(rlam,[]);
    draw.Callback(draw,[]);
    valid = getappdata(fig,'edgepreserve_preview');
    assert(~isempty(valid) && valid.params.rlam==6 && valid.params.adapt_sigma==3);
    assert(valid.report_smooth.is_valid && valid.report_blend.is_valid);
    display = findobj(fig,'Tag','display'); display.Value = 2; display.Callback(display,[]);
    exportapp(fig,fullfile(outputDir,'multiscale_shortest_rlam6_valid_preview.png'));
    rlam.String = '2'; rlam.Callback(rlam,[]);
    draw.Callback(draw,[]);
    retained = getappdata(fig,'edgepreserve_preview');
    assert(isequal(retained.S_smooth.r,valid.S_smooth.r) && retained.params.rlam==6);
    status = strjoin(cellstr(findobj(fig,'Tag','status').String),newline);
    assert(contains(status,'Previous successful preview retained') && contains(status,'PENDING'));
    saveButton = findobj(fig,'Tag','save'); saveButton.Callback(saveButton,[]);
    folders = dir(fullfile(outputDir,source.geometry_name,'preview_*'));
    assert(isscalar(folders));
    directory = fullfile(folders.folder,folders.name);
    stored = load(fullfile(directory,'parameters.mat'));
    assert(stored.parameters.settings.rlam==6 && stored.parameters.settings.adapt_sigma==3);
    manifest = readtable(fullfile(directory,'batch_manifest.csv'),'TextType','string');
    assert(height(manifest)==2 && all(manifest.status=="saved"));
    for k=1:height(manifest)
        saved = surfsmooth3d.surfer.load_from_file(fullfile(directory,manifest.filename(k)));
        assert(saved.npatches==4438 && all(saved.norders==6));
    end
    fprintf('PASS: rlam6 recovery, validity, failed-preview retention and safe pair export with pending rlam2.\n');
end
for k=1:numel(source.source_paths)
    assert(strcmp(surfsmooth3d.edgepreserve.file_sha256(source.source_paths{k}),source.source_hashes{k}));
end
fprintf('PASS: actual shortest-side GUI stop opens full-step diagnostics; source unchanged and no failed save.\n');
end

function close_figures(fig)
if isgraphics(fig), close(fig); end
delete(findall(groot,'Type','figure','Name','Two-stage Newton failure'));
end

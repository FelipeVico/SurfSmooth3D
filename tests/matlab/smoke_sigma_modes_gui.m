function smoke_sigma_modes_gui(packageDir, outputParent)
%SMOKE_SIGMA_MODES_GUI Exercise pending modes, drawn saves, failures and batches.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
if nargin < 2, outputParent = tempname; end
if ~isfolder(outputParent), mkdir(outputParent); end
fig = surfsmooth3d.edgepreserve.launch_gui(root,packageDir,'off',outputParent);
cleanup = onCleanup(@() close_if_open(fig));
selector = findobj(fig,'Tag','adapt');
assert(strcmp(selector.Style,'popupmenu') && selector.Value==3);
assert(isequal(selector.String, {surfsmooth3d.edgepreserve.sigma_mode_label(0); ...
    surfsmooth3d.edgepreserve.sigma_mode_label(3);surfsmooth3d.edgepreserve.sigma_mode_label(1)}));
set(findobj(fig,'Tag','previewOrder'),'String','4');
set(findobj(fig,'Tag','previewLevel'),'String','0');
set(findobj(fig,'Tag','rlam'),'String','4');
draw = findobj(fig,'Tag','draw');
saveButton = findobj(fig,'Tag','save');
status = findall(fig,'Tag','status');
assert(isscalar(status));
source = getappdata(fig,'edgepreserve_source');
draw.Callback(draw,[]);
preview = getappdata(fig,'edgepreserve_preview');
assert(~isempty(preview) && preview.params.adapt_sigma==1);
for entry = [2 1 3]
    previous = getappdata(fig,'edgepreserve_preview');
    selector.Value = entry; selector.Callback(selector,[]);
    current = getappdata(fig,'edgepreserve_preview');
    assert(isequal(current.S_smooth.r,previous.S_smooth.r));
    assert(current.params.adapt_sigma==previous.params.adapt_sigma);
    assert(contains(status_text(status),'PENDING'));
    assert(contains(status_text(status),surfsmooth3d.edgepreserve.sigma_mode_label(previous.params.adapt_sigma)));
    % Saving with a pending mode must keep the drawn surface's provenance.
    saveButton.Callback(saveButton,[]);
    assert_latest_preview(outputParent,source,previous.params.adapt_sigma);
    draw.Callback(draw,[]);
    current = getappdata(fig,'edgepreserve_preview');
    flags = [0 3 1];
    assert(current.params.adapt_sigma==flags(entry));
    assert(contains(status_text(status),'matches expensive controls'));
    if flags(entry)==0
        assert(max(current.sigma)-min(current.sigma)<1e-14);
        assert(all(current.gradSigma==0,'all'));
    end
end
% The selector must reach both full and partial batch computations.
selector.Value = 2; selector.Callback(selector,[]);
set(findobj(fig,'Tag','orders'),'String','[4]');
set(findobj(fig,'Tag','levels'),'String','[0 1]');
batch = findobj(fig,'Tag','batch'); batch.Callback(batch,[]);
folders = dir(fullfile(outputParent,source.geometry_name,'batch_*'));
assert(isscalar(folders));
directory = fullfile(folders(1).folder,folders(1).name);
stored = load(fullfile(directory,'parameters.mat'));
assert(stored.parameters.settings.adapt_sigma==3);
assert(strcmp(stored.parameters.sigma_mode,surfsmooth3d.edgepreserve.sigma_mode_label(3)));
manifest = readtable(fullfile(directory,'batch_manifest.csv'),'TextType','string');
assert(height(manifest)==4 && all(manifest.status=="saved"), ...
    strjoin(string(manifest.error_message),newline));
assert(all(contains(manifest.filename,'_adapt3')));
for k = 1:height(manifest)
    saved = surfsmooth3d.surfer.load_from_file(fullfile(directory,manifest.filename(k)));
    assert(saved.npatches==source.coarse_triangle_count*4^manifest.refinement(k));
    assert(all(saved.norders==4));
end
exportapp(fig,fullfile(outputParent,'sigma_selector_gui.png'));

% Isolate a controlled solver exception without changing production APIs.
mockDirectory = tempname; mkdir(mockDirectory);
mkdir(fullfile(mockDirectory,'+surfsmooth3d'));
fid = fopen(fullfile(mockDirectory,'+surfsmooth3d','multiscale_mesher.m'),'w');
fprintf(fid,['function surfaces = multiscale_mesher(varargin)\n' ...
    'error(''sigma_modes:ControlledFailure'',''Deliberate recompute failure.'');\nend\n']);
fclose(fid);
addpath(mockDirectory,'-begin'); clear surfsmooth3d.multiscale_mesher;
restore = onCleanup(@() remove_mock(mockDirectory));
previous = getappdata(fig,'edgepreserve_preview');
draw.Callback(draw,[]);
retained = getappdata(fig,'edgepreserve_preview');
assert(isequal(retained.S_smooth.r,previous.S_smooth.r));
assert(retained.params.adapt_sigma==previous.params.adapt_sigma);
assert(contains(status_text(status),'Deliberate recompute failure'));
assert(contains(status_text(status),'PENDING'));
for k = 1:numel(source.source_paths)
    assert(strcmp(surfsmooth3d.edgepreserve.file_sha256(source.source_paths{k}),source.source_hashes{k}));
end
fprintf('PASS: sigma selector, pending/drawn modes, saves, shortest-side batch and failed-preview retention.\n');
end

function assert_latest_preview(outputParent,source,flag)
folders = dir(fullfile(outputParent,source.geometry_name,'preview_*'));
[~,index] = max([folders.datenum]);
directory = fullfile(folders(index).folder,folders(index).name);
stored = load(fullfile(directory,'parameters.mat'));
assert(stored.parameters.settings.adapt_sigma==flag);
assert(strcmp(stored.parameters.sigma_mode,surfsmooth3d.edgepreserve.sigma_mode_label(flag)));
manifest = readtable(fullfile(directory,'batch_manifest.csv'),'TextType','string');
assert(height(manifest)==2 && all(manifest.status=="saved"));
assert(all(contains(manifest.filename,sprintf('_adapt%d',flag))));
end

function close_if_open(fig)
if isgraphics(fig), close(fig); end
end

function text = status_text(control)
text = strjoin(cellstr(control.String),newline);
end

function remove_mock(directory)
rmpath(directory); clear surfsmooth3d.multiscale_mesher;
rmdir(directory,'s');
end

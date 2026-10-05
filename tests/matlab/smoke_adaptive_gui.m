function smoke_adaptive_gui(packageDir,outputParent)
%SMOKE_ADAPTIVE_GUI Real CAD case, all view modes, pending edits and save.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
if nargin < 1 || isempty(packageDir)
    packageDir = fullfile(fullfile(root,'tests','fixtures','cad'), ...
        'convergence_intersection_two_balls','p4_mf0p10000_rl0');
end
if nargin < 2, outputParent = tempname; mkdir(outputParent); end
if ~isfolder(outputParent), mkdir(outputParent); end
fig = surfsmooth3d.adaptivesmoother.launch_gui(root,packageDir,'off',outputParent);
cleanup = onCleanup(@() delete(fig));
assert(isempty(getappdata(fig,'adaptive_result')));
set(findobj(fig,'Tag','order'),'String','4');
set(findobj(fig,'Tag','tolerance'),'String','1e-3');
set(findobj(fig,'Tag','depth'),'String','2');
invoke('draw');
result = getappdata(fig,'adaptive_result');
assert(~isempty(result));
save(fullfile(outputParent,'adaptive_real_result.mat'),'result');
pkg = getappdata(fig,'adaptive_source');
assert(result.surface.npatches <= pkg.coarse_triangle_count*4^2);
ax = findobj(fig,'Type','axes'); ax = ax(1);
cb = findobj(fig,'Type','colorbar');
limits = [ax.XLim;ax.YLim;ax.ZLim]; position = cb.Position;
for mode = [3 2 1]
    set(findobj(fig,'Tag','display'),'Value',mode); invoke('display');
    for color=1:7
        set(findobj(fig,'Tag','color'),'Value',color); invoke('color');
        assert(isequal(limits,[ax.XLim;ax.YLim;ax.ZLim]) && isequal(cb.Position,position));
    end
end
set(findobj(fig,'Tag','color'),'Value',2); invoke('color');
set(findobj(fig,'Tag','wire'),'Value',1); invoke('wire');
exportapp(fig,fullfile(outputParent,'adaptive_gui.png'));
set(findobj(fig,'Tag','rlam'),'String','3'); invoke('rlam');
same = getappdata(fig,'adaptive_result');
assert(isequal(same.surface.r,result.surface.r) && same.settings.rlam==2);
status = findobj(fig,'Tag','status');
assert(any(contains(string(status.String),'Pending')));
if result.report.is_valid
    invoke('save'); filename = getappdata(fig,'adaptive_saved_file');
    assert(isfile(filename));
    stored = load(fullfile(fileparts(filename),'parameters.mat'));
    assert(stored.parameters.drawn_settings.rlam==2);
end
% Deliberately over-smooth: the level set disappears and Newton must fail safely.
set(findobj(fig,'Tag','rlam'),'String','0.05'); invoke('draw');
same = getappdata(fig,'adaptive_result');
assert(isequal(same.surface.r,result.surface.r));
assert(any(contains(string(status.String),'Previous successful preview retained.')));
set(findobj(fig,'Tag','order'),'String','99'); invoke('draw');
same = getappdata(fig,'adaptive_result');
assert(isequal(same.surface.r,result.surface.r));
for k=1:numel(pkg.source_paths)
    assert(strcmp(surfsmooth3d.edgepreserve.file_sha256(pkg.source_paths{k}),pkg.source_hashes{k}));
end
save(fullfile(outputParent,'adaptive_real_result.mat'),'result');
fprintf('PASS: adaptive GUI views, fixed layout, pending settings, prior-preview preservation and source hashes.\n');
fprintf('Real case: %d patches (uniform at deepest level: %d), %d unresolved; validity %d.\n', ...
    result.surface.npatches,pkg.coarse_triangle_count*4^result.info.achieved_depth, ...
    result.info.unresolved_count,result.report.is_valid);

    function invoke(tag)
        h = findobj(fig,'Tag',tag); callback = h.Callback; callback(h,[]);
    end
end

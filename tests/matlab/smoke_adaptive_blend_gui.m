function smoke_adaptive_blend_gui(packageDir,output)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
lockBaseline = adaptive_blend_lock_baseline();
if nargin<1 || isempty(packageDir)
    packageDir = fullfile(fullfile(root,'tests','fixtures','cad'), ...
        'convergence_intersection_two_balls','p4_mf0p10000_rl0');
end
if nargin<2, output = fullfile(tempdir,'adaptive_blend_gui'); end
if ~isfolder(output), mkdir(output); end
f = surfsmooth3d.adaptiveblend.launch_gui(root,packageDir,'off',output,500000);
cleanup = onCleanup(@() close(f));
pkg = getappdata(f,'adaptiveblend_source'); assert(~isempty(pkg));
assert(isempty(getappdata(f,'adaptiveblend_result')));
set_value('order','4'); set_value('tolerance','3e-3'); set_value('depth','2');
set_value('edges',mat2str([pkg.edges.edge_id]));
invoke('draw'); result = getappdata(f,'adaptiveblend_result'); assert(~isempty(result));
assert(result.report.is_valid);
save(fullfile(output,'result.mat'),'result','pkg');
surfsmooth3d.edgepreserve.write_point_source(pkg.surface,fullfile(output,'cad_points.txt'));
copyfile(pkg.scaffold_file,fullfile(output,'scaffold.gidmsh'));
ax = findobj(f,'Type','axes'); cb = findobj(f,'Type','colorbar');
limits = [ax.XLim ax.YLim ax.ZLim]; axpos = ax.Position; cbpos = cb.Position;
for display = 1:4
    h = findobj(f,'Tag','display'); h.Value = display;
    for color = 1:8
        h = findobj(f,'Tag','color'); h.Value = color; invoke('color');
        assert(max(abs([ax.XLim ax.YLim ax.ZLim]-limits))<1e-10);
        assert(max(abs(ax.Position-axpos))<1e-10 && max(abs(cb.Position-cbpos))<1e-10);
    end
end
h = findobj(f,'Tag','display'); h.Value = 1;
h = findobj(f,'Tag','color'); h.Value = 2;
h = findobj(f,'Tag','labels'); h.Value = 0;
h = findobj(f,'Tag','showedges'); h.Value = 0;
h = findobj(f,'Tag','wire'); h.Value = 1; invoke('wire');
exportapp(f,fullfile(output,'adaptive_blend_gui.png'));
set_value('crad','3.1'); invoke('crad');
assert(any(contains(string(findobj(f,'Tag','status').String),'Pending settings')));
invoke('save'); manifest = getappdata(f,'adaptiveblend_saved'); assert(~isempty(manifest));
saved = load(fullfile(fileparts(manifest.filename(1)),'parameters.mat'));
assert(saved.parameters.drawn_settings.c_rad==2.7);
assert(isequal(getappdata(f,'adaptiveblend_result').surface.r,result.surface.r));
set_value('order','99'); invoke('draw');
assert(isequal(getappdata(f,'adaptiveblend_result').surface.r,result.surface.r));
set_value('order','4'); set_value('rlam','0.05'); invoke('draw');
assert(isequal(getappdata(f,'adaptiveblend_result').surface.r,result.surface.r));
assert(mislocked('surfsmooth3d_adaptive_blend_routs')==lockBaseline);
for j = 1:numel(pkg.source_paths)
    assert(strcmp(pkg.source_hashes{j},surfsmooth3d.edgepreserve.file_sha256(pkg.source_paths{j})));
end
fprintf('PASS: real adaptive blend, GUI modes, pending state, guarded failure and saved provenance.\n');

    function set_value(tag,value)
        control = findobj(f,'Tag',tag); control.String = value;
    end
    function invoke(tag)
        h = findobj(f,'Tag',tag); callback = h.Callback; callback(h,[]); drawnow;
    end
end

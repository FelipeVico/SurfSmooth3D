function smoke_edgepreserve_real_package(packageDir, outputParent)
%SMOKE_EDGEPRESERVE_REAL_PACKAGE Real Fortran solve, fixed sigma and GUI smoke.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
if nargin < 1 || isempty(packageDir)
    packageDir = fullfile(fullfile(root,'tests','fixtures','cad'), ...
        'convergence_intersection_two_balls','p8_mf0p10000_rl1');
end
if nargin < 2, outputParent = tempname; mkdir(outputParent); end
pkg = surfsmooth3d.edgepreserve.load_package(packageDir);
[work,cleanup] = surfsmooth3d.edgepreserve.stage_source(pkg); %#ok<ASGLU>
settings = surfsmooth3d.edgepreserve.default_settings();
settings.orders = [4 6]; settings.levels = [0 1];
settings.selectedEdgeIds = pkg.edges(1).edge_id;
opts = surfsmooth3d.edgepreserve.mesher_options(pkg,work,settings,0);
targets = pkg.surface.r(:,1:30);
for order = [4 6 8]
    [sigma,gradient] = surfsmooth3d.multiscale_mesher_sigma_eval(work.scaffold_file,order,opts,targets);
    if order == 4
        sigma0 = sigma; gradient0 = gradient;
    else
        assert(max(abs(sigma-sigma0))<1e-13);
        assert(max(abs(gradient(:)-gradient0(:)))<1e-13);
    end
end
[manifest,directory] = surfsmooth3d.edgepreserve.run_batch(pkg,settings,outputParent);
assert(height(manifest)==8);
for k = 1:height(manifest)
    if manifest.status(k) ~= "saved", continue; end
    S = surfsmooth3d.surfer.load_from_file(fullfile(directory,manifest.filename(k)));
    assert(all(S.norders==manifest.order(k)));
    assert(S.npatches==pkg.coarse_triangle_count*4^manifest.refinement(k));
end
assert(all(manifest.status=="saved"),strjoin(manifest.error_message,newline));
for k = 1:numel(pkg.source_paths)
    assert(strcmp(surfsmooth3d.edgepreserve.file_sha256(pkg.source_paths{k}),pkg.source_hashes{k}));
end
fig = surfsmooth3d.edgepreserve.launch_gui(root,packageDir,'off',outputParent);
closer = onCleanup(@() delete_if_open(fig));
assert(isempty(getappdata(fig,'edgepreserve_preview')));
set(findobj(fig,'Tag','previewOrder'),'String','4');
set(findobj(fig,'Tag','previewLevel'),'String','0');
invoke(fig,'draw');
result = getappdata(fig,'edgepreserve_preview');
assert(~isempty(result) && result.order==4 && result.level==0);
axesHandle = findobj(fig,'Type','axes');
axesHandle = axesHandle(1);
limits = [axesHandle.XLim;axesHandle.YLim;axesHandle.ZLim];
bar = findobj(fig,'Type','colorbar');
barPosition = bar.Position;
for mode = [3 2 1 4]
    set(findobj(fig,'Tag','display'),'Value',mode); invoke(fig,'display');
    assert(isequal(limits,[axesHandle.XLim;axesHandle.YLim;axesHandle.ZLim]));
    assert(isequal(barPosition,bar.Position));
end
set(findobj(fig,'Tag','display'),'Value',2); invoke(fig,'display');
set(findobj(fig,'Tag','color'),'Value',2); invoke(fig,'color');
assert(axesHandle.CLim(1)<0 && axesHandle.CLim(2)>0);
set(findobj(fig,'Tag','rlam'),'String','3'); invoke(fig,'rlam');
unchanged = getappdata(fig,'edgepreserve_preview');
assert(isequal(unchanged.S_smooth.r,result.S_smooth.r) && unchanged.params.rlam==2);
invoke(fig,'save');
set(findobj(fig,'Tag','display'),'Value',1); invoke(fig,'display');
set(findobj(fig,'Tag','color'),'Value',4); invoke(fig,'color');
exportapp(fig,fullfile(outputParent,'gui_smoke.png'));
close(fig);
fprintf('PASS: real package, sigma invariance, eight exports, GUI switching/pending settings.\n%s\n',directory);
end

function invoke(fig,tag)
h = findobj(fig,'Tag',tag);
callback = h.Callback;
callback(h,[]);
end

function delete_if_open(fig)
if isgraphics(fig), close(fig); end
end

function fig = launch_gui(repoRoot, packageDir, visibility, outputParent, maxPoints)
%LAUNCH_GUI Independent fully smooth adaptive workflow; no blend controls.
if nargin < 2, packageDir = ''; end
if nargin < 3, visibility = 'on'; end
if nargin < 4 || isempty(outputParent)
    outputParent = fullfile(repoRoot,'outputs','smoother');
end
if nargin < 5, maxPoints = 2000000; end
sigmaFlags = [0 3 1];
state = struct('pkg',[],'result',[],'busy',false,'limits',[], ...
    'notice','Choose a geometry package.','failureFigure',[]);
inputRoot = surfsmooth3d.edgepreserve.default_input_root(repoRoot);
fig = figure('Name','Two-stage adaptive surface smoother','NumberTitle','off', ...
    'Position',[90 65 1420 920],'Visible',visibility,'CloseRequestFcn',@close_viewer);
ax = axes('Parent',fig,'Units','normalized','Position',[.05 .19 .53 .71]);
cb = colorbar(ax);
set(cb,'Units','normalized','Position',[.605 .19 .017 .71]);
panel = uipanel(fig,'Title','Fully smooth / adaptive refinement','Units','normalized', ...
    'Position',[.665 .025 .325 .95]);
status = uicontrol(fig,'Style','text','Units','normalized','Tag','status', ...
    'Position',[.035 .015 .595 .145],'HorizontalAlignment','left','FontSize',11);
controls = struct();
control('pushbutton','Browse Package',[.04 .935 .44 .045],@browse,'browse');
control('pushbutton','Output Folder',[.52 .935 .44 .045],@output,'output');
controls.source = control('text','No package loaded.',[.04 .85 .92 .075],[],'source');
field('Output order','6',.78,'order');
field('rlam','2',.70,'rlam');
control('text','Sigma mode',[.04 .655 .92 .025],[],'');
controls.sigma = control('popupmenu',arrayfun(@surfsmooth3d.edgepreserve.sigma_mode_label, ...
    sigmaFlags,'UniformOutput',false),[.04 .615 .92 .037],@pending,'sigma');
controls.sigma.Value = 3;
field('Adaptive tolerance','1e-4',.535,'tolerance');
field('Maximum refinement depth','3',.455,'depth');
controls.draw = control('pushbutton','Draw / Recompute',[.04 .386 .92 .05],@draw,'draw');
controls.save = control('pushbutton','Save Smooth Surface',[.04 .321 .92 .05],@save_surface,'save');
control('text','Display mode',[.04 .280 .92 .025],[],'');
controls.display = control('popupmenu',{'fully smooth','CAD source','original scaffold'}, ...
    [.04 .242 .92 .035],@view_changed,'display');
control('text','Color mode',[.04 .200 .92 .025],[],'');
controls.color = control('popupmenu',{'patch resolution indicator','refinement depth', ...
    'unresolved patches','signed mean curvature','absolute mean curvature', ...
    'jacobian ratio to CAD','geometry health'},[.04 .162 .92 .035],@view_changed,'color');
controls.wire = control('checkbox','Show patch boundaries',[.04 .098 .92 .035],@view_changed,'wire');
rotation = rotate3d(fig);
rotation.Enable = 'on';
rotation.ActionPostCallback = @(~,~) light_update();
if ~isempty(packageDir)
    load_source(packageDir);
elseif strcmp(visibility,'on')
    browse();
end
update_status();

    function h = control(style,text,position,callback,tag)
        h = uicontrol(panel,'Style',style,'String',text,'Units','normalized', ...
            'Position',position,'Tag',tag,'HorizontalAlignment','left');
        if ~isempty(callback), h.Callback = callback; end
    end
    function field(label,value,y,tag)
        control('text',label,[.04 y+.037 .92 .025],[],'');
        controls.(tag) = control('edit',value,[.04 y .92 .037],@pending,tag);
    end
    function s = read_settings()
        s = struct('order',str2double(controls.order.String),'rlam',str2double(controls.rlam.String), ...
            'adapt_sigma',sigmaFlags(controls.sigma.Value),'eps_adapt',str2double(controls.tolerance.String), ...
            'max_refine',str2double(controls.depth.String),'max_points',maxPoints);
        validateattributes(s.order,{'double'},{'scalar','integer','>=',1,'<=',20});
        validateattributes(s.rlam,{'double'},{'scalar','finite','positive'});
        validateattributes(s.eps_adapt,{'double'},{'scalar','finite','>',0,'<',1});
        validateattributes(s.max_refine,{'double'},{'scalar','integer','>=',0,'<=',20});
    end
    function browse(varargin)
        if state.busy, return; end
        start = inputRoot;
        if ~isempty(state.pkg), start = state.pkg.package_dir; end
        folder = uigetdir(start,'Choose geometry package');
        if ~isequal(folder,0), load_source(folder); end
    end
    function load_source(folder)
        try
            pkg = surfsmooth3d.edgepreserve.load_package(folder);
            state.pkg = pkg;
            state.result = [];
            state.limits = [];
            expand_limits([pkg.scaffold.nodes pkg.surface.r]);
            controls.source.String = sprintf('%s / %s\nCAD order %d, source refinement %d\nCoarse triangles %d', ...
                pkg.geometry_name,pkg.name,pkg.source_order,pkg.source_level,pkg.coarse_triangle_count);
            state.notice = 'Fixed CAD source loaded. Preview not computed.';
            setappdata(fig,'adaptive_source',pkg);
            setappdata(fig,'adaptive_result',[]);
            refresh_plot();
        catch exception
            show_error(exception);
        end
    end
    function output(varargin)
        if state.busy, return; end
        folder = uigetdir(outputParent,'Choose output parent folder');
        if ~isequal(folder,0), outputParent = folder; end
    end
    function pending(varargin)
        if state.busy, return; end
        update_status();
    end
    function draw(varargin)
        if state.busy || isempty(state.pkg), return; end
        state.busy = true;
        cleanup = onCleanup(@release);
        try
            settings = read_settings();
            state.notice = 'Running adaptive two-stage smoother...';
            update_status(); drawnow;
            result = surfsmooth3d.adaptivesmoother.solve(state.pkg,settings);
            state.result = result;
            state.notice = result.info.stop_reason;
            expand_limits(result.surface.r);
            setappdata(fig,'adaptive_result',result);
            refresh_plot();
        catch exception
            show_error(exception);
        end
    end
    function release()
        state.busy = false; %#ok<MOCUP> Shared workspace retained by figure callbacks.
        if isgraphics(fig), update_status(); end
    end
    function save_surface(varargin)
        if state.busy || isempty(state.result), return; end
        state.busy = true;
        cleanup = onCleanup(@release);
        try
            filename = surfsmooth3d.adaptivesmoother.save_result(state.pkg,state.result,outputParent);
            state.notice = ['Saved drawn settings: ' filename];
            setappdata(fig,'adaptive_saved_file',filename);
        catch exception
            show_error(exception);
        end
    end
    function show_error(exception)
        [report,filename] = surfsmooth3d.edgepreserve.capture_newton_failure(exception);
        if ~isempty(report)
            fprintf('Newton diagnostic preserved: %s\n',filename);
            state.failureFigure = surfsmooth3d.edgepreserve.plot_newton_failure( ...
                state.pkg,report,state.failureFigure,visibility);
        end
        state.notice = ['Error: ' exception.message];
        if ~isempty(state.result)
            state.notice = [state.notice ' Previous successful preview retained.'];
        end
        warning('%s',exception.message);
        update_status();
    end
    function view_changed(varargin)
        if ~state.busy, refresh_plot(); end
    end
    function expand_limits(r)
        lower = min(r,[],2); upper = max(r,[],2);
        padding = .04*max(max(upper-lower),eps);
        limits = [lower-padding upper+padding];
        if isempty(state.limits)
            state.limits = limits;
        else
            state.limits = [min(state.limits(:,1),limits(:,1)) max(state.limits(:,2),limits(:,2))];
        end
    end
    function refresh_plot()
        if isempty(state.pkg), return; end
        [azimuth,elevation] = view(ax);
        angles = [azimuth elevation];
        cla(ax); hold(ax,'on');
        label = 'constant';
        values = [0 1];
        signed = false;
        if controls.display.Value == 3
            surfsmooth3d.edgepreserve.plot_gidmsh(ax,state.pkg.scaffold,'solid');
            titleText = 'original scaffold';
        else
            S = state.pkg.surface;
            result = [];
            titleText = 'CAD source';
            if controls.display.Value == 1 && ~isempty(state.result)
                result = state.result; S = result.surface;
                titleText = 'fully smooth (drawn settings)';
            end
            values = zeros(S.npts,1);
            switch controls.color.Value
                case 1
                    if ~isempty(result)
                        values = log10(max(max(result.info.indicators,[],2),1e-16));
                        label = 'log10 maximum resolution indicator';
                    end
                case 2
                    if ~isempty(result), values = result.info.depth; label = 'refinement depth'; end
                case 3
                    if ~isempty(result), values = double(result.info.unresolved); label = 'unresolved'; end
                case 4
                    values = S.mean_curv; label = 'signed mean curvature'; signed = true;
                case 5
                    values = abs(S.mean_curv); label = 'absolute mean curvature';
                case 6
                    if ~isempty(result), values = result.report.jac_ratio_to_ref; label = 'Jacobian / CAD'; end
                case 7
                    if ~isempty(result), values = result.report.failure_code; label = 'geometry failure code'; end
            end
            axes(ax);
            surfsmooth3d.edgepreserve.plot_surface(S,values,{'EdgeColor','none'},max(6,double(S.norders(1))+2));
            if controls.wire.Value, patch_boundaries(S); end
        end
        finite = values(isfinite(values));
        if isempty(finite), limits = [0 1]; else, limits = [min(finite(:)) max(finite(:))]; end
        if signed, limits = [-1 1]*max(max(abs(limits)),eps); end
        if limits(1)==limits(2), limits = limits+[-.5 .5]*max(1,abs(limits(1))); end
        clim(ax,limits); cb.Label.String = label;
        axis(ax,'equal'); xlim(ax,state.limits(1,:)); ylim(ax,state.limits(2,:)); zlim(ax,state.limits(3,:));
        if isequal(angles,[0 90]), view(ax,3); else, view(ax,angles); end
        title(ax,{[state.pkg.geometry_name ': ' titleText],label},'Interpreter','none');
        xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
        ax.Position = [.05 .19 .53 .71]; cb.Position = [.605 .19 .017 .71];
        light_update(); update_status();
    end
    function patch_boundaries(S)
        t = linspace(0,1,20);
        uv = [t 1-t zeros(size(t)); zeros(size(t)) t 1-t];
        pol = surfsmooth3d.internal.koorn.pols(double(S.norders(1)),uv);
        xyz = nan(3,size(uv,2)+1,S.npatches);
        for k = 1:S.npatches, xyz(:,1:end-1,k) = S.srccoefs{k}(1:3,:)*pol; end
        xyz = reshape(xyz,3,[]);
        plot3(ax,xyz(1,:),xyz(2,:),xyz(3,:),'Color',[.2 .2 .2],'LineWidth',.4);
    end
    function light_update()
        delete(findobj(ax,'Type','light'));
        camlight(ax,'headlight'); lighting(ax,'gouraud');
    end
    function update_status()
        if ~isgraphics(fig), return; end
        text = state.notice;
        if ~isempty(state.result)
            result = state.result; info = result.info;
            health = 'PASS';
            if result.report.has_warnings, health = 'WARN'; end
            if ~result.report.is_valid, health = 'FAIL'; end
            text = sprintf('%s\nDrawn: p=%d, rlam=%.5g, %s\n%d patches, %d nodes; depth %d; unresolved %d; validity %s', ...
                text,result.settings.order,result.settings.rlam, ...
                surfsmooth3d.edgepreserve.sigma_mode_label(result.settings.adapt_sigma), ...
                result.surface.npatches,result.surface.npts,info.achieved_depth,info.unresolved_count,health);
            text = sprintf('%s\nMax indicators [position tail, normal tail, position check, normal check]: %s', ...
                text,mat2str(max(info.indicators,[],1),3));
            try
                different = ~isequal(read_settings(),result.settings);
            catch
                different = true;
            end
            if different, text = [text newline 'Pending settings not drawn. Save uses the last computed smooth surface.']; end
        end
        status.String = text;
        controls.draw.Enable = on_off(~state.busy && ~isempty(state.pkg));
        controls.save.Enable = on_off(~state.busy && ~isempty(state.result));
    end
    function close_viewer(varargin)
        if state.busy, return; end
        delete(fig);
    end
end

function value = on_off(condition)
value = 'off';
if condition, value = 'on'; end
end

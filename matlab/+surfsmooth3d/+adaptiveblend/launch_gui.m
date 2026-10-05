function fig = launch_gui(repoRoot,packageDir,visibility,outputParent,maxPoints)
%LAUNCH_GUI Adaptive refinement driven by the partially smoothed surface.
if nargin<2, packageDir = ''; end
if nargin<3, visibility = 'on'; end
if nargin<4 || isempty(outputParent)
    outputParent = fullfile(repoRoot,'outputs','smoother');
end
if nargin<5, maxPoints = 2000000; end
sigmaFlags = [0 3 1];
state = struct('pkg',[],'result',[],'busy',false,'limits',[], ...
    'notice','Choose a geometry package.','failureFigure',[]);
inputRoot = surfsmooth3d.edgepreserve.default_input_root(repoRoot);
fig = figure('Name','Adaptive partially smoothed surface','NumberTitle','off', ...
    'Position',[90 65 1420 940],'Visible',visibility,'CloseRequestFcn',@close_viewer);
ax = axes(fig,'Units','normalized','Position',[.05 .20 .53 .70]);
cb = colorbar(ax); cb.Units = 'normalized'; cb.Position = [.605 .20 .017 .70];
panel = uipanel(fig,'Title','Partially smoothed / adaptive refinement', ...
    'Units','normalized','Position',[.665 .025 .325 .95]);
status = uicontrol(fig,'Style','text','Units','normalized','Tag','status', ...
    'Position',[.035 .015 .595 .16],'HorizontalAlignment','left','FontSize',10);
controls = struct();
control('pushbutton','Browse Package',[.04 .935 .44 .045],@browse,'browse');
control('pushbutton','Output Folder',[.52 .935 .44 .045],@output,'output');
controls.source = control('text','No package loaded.',[.04 .85 .92 .075],[],'source');
field('Smooth edge IDs','[]',.79,'edges',.04,.92);
controls.labels = control('checkbox','Show CAD edge labels',[.04 .747 .92 .032],@refresh,'labels');
controls.labels.Value = 1;
field('Output order','6',.687,'order',.04,.43);
field('rlam','2',.687,'rlam',.53,.43);
control('text','Sigma mode',[.04 .648 .92 .025],[],'');
controls.sigma = control('popupmenu',arrayfun(@surfsmooth3d.edgepreserve.sigma_mode_label, ...
    sigmaFlags,'UniformOutput',false),[.04 .609 .92 .037],@pending,'sigma');
controls.sigma.Value = 3;
field('c_eps','0.5',.538,'ceps',.04,.27);
field('c_rad','2.7',.538,'crad',.365,.27);
field('c_tau','1.05',.538,'ctau',.69,.27);
field('Adaptive tolerance','1e-4',.463,'tolerance',.04,.43);
field('Maximum depth','3',.463,'depth',.53,.43);
controls.draw = control('pushbutton','Draw / Recompute',[.04 .391 .92 .05],@draw,'draw');
controls.save = control('pushbutton','Save Current Pair',[.04 .326 .92 .05],@save_pair,'save');
control('text','Display mode',[.04 .282 .92 .025],[],'');
controls.display = control('popupmenu',{'partially smoothed','fully smooth companion', ...
    'CAD source','original scaffold'},[.04 .245 .92 .035],@refresh,'display');
control('text','Color mode',[.04 .202 .92 .025],[],'');
controls.color = control('popupmenu',{'patch resolution indicator','refinement depth', ...
    'unresolved patches','beta','signed mean curvature','absolute mean curvature', ...
    'jacobian ratio to CAD','geometry health'},[.04 .165 .92 .035],@refresh,'color');
controls.wire = control('checkbox','Show patch boundaries',[.04 .107 .92 .032],@refresh,'wire');
controls.showedges = control('checkbox','Show CAD edges',[.04 .060 .92 .032],@refresh,'showedges');
controls.showedges.Value = 1;
rotation = rotate3d(fig); rotation.Enable = 'on';
rotation.ActionPostCallback = @(~,~) light_update();
if ~isempty(packageDir), load_source(packageDir); elseif strcmp(visibility,'on'), browse(); end
update_status();

    function h = control(style,text,position,callback,tag)
        h = uicontrol(panel,'Style',style,'String',text,'Units','normalized', ...
            'Position',position,'Tag',tag,'HorizontalAlignment','left');
        if ~isempty(callback), h.Callback = callback; end
    end
    function field(label,value,y,tag,x,width)
        control('text',label,[x y+.037 width .025],[],'');
        controls.(tag) = control('edit',value,[x y width .037],@pending,tag);
    end
    function s = read_settings()
        s = surfsmooth3d.adaptiveblend.default_settings();
        s.order = str2double(controls.order.String); s.rlam = str2double(controls.rlam.String);
        s.adapt_sigma = sigmaFlags(controls.sigma.Value);
        s.eps_adapt = str2double(controls.tolerance.String);
        s.max_refine = str2double(controls.depth.String); s.max_points = maxPoints;
        s.c_eps = str2double(controls.ceps.String); s.c_rad = str2double(controls.crad.String);
        s.c_tau = str2double(controls.ctau.String);
        s.selectedEdgeIds = surfsmooth3d.edgepreserve.parse_integer_list(controls.edges.String,1,flintmax,true);
        surfsmooth3d.adaptiveblend.validate_settings(s,state.pkg);
    end
    function browse(varargin)
        if state.busy, return; end
        start = inputRoot; if ~isempty(state.pkg), start = state.pkg.package_dir; end
        folder = uigetdir(start,'Choose geometry package');
        if ~isequal(folder,0), load_source(folder); end
    end
    function load_source(folder)
        try
            pkg = surfsmooth3d.edgepreserve.load_package(folder);
            state.pkg = pkg; state.result = []; state.limits = [];
            controls.edges.String = '[]';
            expand_limits([pkg.surface.r pkg.scaffold.nodes]);
            controls.source.String = sprintf('%s / %s\nCAD order %d, source refinement %d\nCoarse triangles %d, CAD edges %d', ...
                pkg.geometry_name,pkg.name,pkg.source_order,pkg.source_level,pkg.coarse_triangle_count,numel(pkg.edges));
            state.notice = 'Fixed CAD loaded. Empty edge selection gives a CAD-only partial surface.';
            setappdata(fig,'adaptiveblend_source',pkg); setappdata(fig,'adaptiveblend_result',[]);
            refresh();
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
        if ~state.busy, update_status(); end
    end
    function draw(varargin)
        if state.busy || isempty(state.pkg), return; end
        state.busy = true; cleanup = onCleanup(@release);
        try
            s = read_settings();
            state.notice = 'Computing adaptive partial surface...'; update_status(); drawnow;
            result = surfsmooth3d.adaptiveblend.solve(state.pkg,s);
            state.result = result; state.notice = result.info.stop_reason;
            expand_limits([result.surface.r result.smooth.r]);
            setappdata(fig,'adaptiveblend_result',result);
            refresh_plot();
        catch exception
            show_error(exception);
        end
    end
    function release()
        state.busy = false; %#ok<MOCUP> Figure callbacks share this state.
        if isgraphics(fig), update_status(); end
    end
    function save_pair(varargin)
        if state.busy || isempty(state.result), return; end
        state.busy = true; cleanup = onCleanup(@release);
        try
            manifest = surfsmooth3d.adaptiveblend.save_result(state.pkg,state.result,outputParent);
            state.notice = sprintf('Saved %d/2 surfaces using drawn settings.',sum(manifest.status=="saved"));
            setappdata(fig,'adaptiveblend_saved',manifest);
        catch exception
            show_error(exception);
        end
    end
    function show_error(exception)
        [report,filename] = surfsmooth3d.edgepreserve.capture_newton_failure(exception);
        if ~isempty(report)
            fprintf('Newton diagnostic preserved: %s\n',filename);
            state.failureFigure = surfsmooth3d.edgepreserve.plot_newton_failure(state.pkg,report,state.failureFigure,visibility);
        end
        state.notice = ['Error: ' exception.message];
        if ~isempty(state.result), state.notice = [state.notice ' Previous successful preview retained.']; end
        warning('%s',exception.message); update_status();
    end
    function expand_limits(r)
        lower = min(r,[],2); upper = max(r,[],2);
        padding = .04*max(max(upper-lower),eps); limits = [lower-padding upper+padding];
        if isempty(state.limits), state.limits = limits;
        else, state.limits = [min(state.limits(:,1),limits(:,1)) max(state.limits(:,2),limits(:,2))]; end
    end
    function refresh(varargin)
        if ~state.busy, refresh_plot(); end
    end
    function refresh_plot()
        if isempty(state.pkg), return; end
        [az,el] = view(ax); cla(ax); hold(ax,'on');
        label = 'constant'; values = [0 1]; signed = false;
        if controls.display.Value==4
            surfsmooth3d.edgepreserve.plot_gidmsh(ax,state.pkg.scaffold,'solid'); titleText = 'original scaffold';
        else
            S = state.pkg.surface; titleText = 'CAD source'; result = []; report = [];
            if controls.display.Value<=2 && ~isempty(state.result)
                result = state.result; S = result.surface; report = result.report;
                titleText = 'partially smoothed (drawn settings)';
                if controls.display.Value==2
                    S = result.smooth; report = result.report_smooth;
                    titleText = 'fully smooth companion (not tolerance-certified)';
                end
            end
            values = zeros(S.npts,1);
            switch controls.color.Value
                case 1
                    if ~isempty(result)
                        values = log10(max(max(result.info.indicators,[],2),1e-16));
                        finite = values(isfinite(values)); cap = 0;
                        if ~isempty(finite), cap = max(0,max(finite)); end
                        values(~isfinite(values)) = cap+1;
                        label = 'log10 partial-patch indicator (invalid checks at top)';
                    end
                case 2
                    if ~isempty(result), values = result.info.depth; label = 'refinement depth'; end
                case 3
                    if ~isempty(result), values = double(result.info.unresolved); label = 'partial patches unresolved'; end
                case 4
                    if ~isempty(result) && controls.display.Value==1, values = result.beta; label = 'smooth weight beta'; end
                case 5
                    values = S.mean_curv; label = 'signed mean curvature'; signed = true;
                case 6
                    values = abs(S.mean_curv); label = 'absolute mean curvature';
                case 7
                    if ~isempty(report), values = report.jac_ratio_to_ref; label = 'Jacobian / CAD'; end
                case 8
                    if ~isempty(report), values = report.failure_code; label = 'geometry failure code'; end
            end
            axes(ax);
            surfsmooth3d.edgepreserve.plot_surface(S,values,{'EdgeColor','none'},max(6,double(S.norders(1))+2));
            if controls.wire.Value, patch_boundaries(S); end
        end
        overlay_edges();
        finite = values(isfinite(values));
        if isempty(finite), limits = [0 1]; else, limits = [min(finite(:)) max(finite(:))]; end
        if signed, limits = [-1 1]*max(max(abs(limits)),eps); end
        if limits(1)==limits(2), limits = limits+[-.5 .5]*max(1,abs(limits(1))); end
        if controls.color.Value==4, limits = [0 1]; end
        clim(ax,limits); cb.Label.String = label;
        axis(ax,'equal'); xlim(ax,state.limits(1,:)); ylim(ax,state.limits(2,:)); zlim(ax,state.limits(3,:));
        if isequal([az el],[0 90]), view(ax,3); else, view(ax,az,el); end
        title(ax,{[state.pkg.geometry_name ': ' titleText],label},'Interpreter','none');
        xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
        ax.Position = [.05 .20 .53 .70]; cb.Position = [.605 .20 .017 .70];
        light_update(); update_status();
    end
    function overlay_edges()
        if ~controls.showedges.Value && ~controls.labels.Value, return; end
        selected = [];
        if ~isempty(state.result), selected = state.result.settings.selectedEdgeIds; end
        for edge = state.pkg.edges(:).'
            xyz = cat(2,edge.panels.xyz);
            if controls.showedges.Value
                width = .7; color = [.4 .4 .4];
                if ismember(edge.edge_id,selected), width = 2; color = [0 0 0]; end
                for piece = edge.panels(:).'
                    plot3(ax,piece.xyz(1,:),piece.xyz(2,:),piece.xyz(3,:),'Color',color,'LineWidth',width);
                end
            end
            if controls.labels.Value && ~isempty(xyz)
                point = xyz(:,ceil(size(xyz,2)/2));
                text(ax,point(1),point(2),point(3),sprintf('E%d',edge.edge_id), ...
                    'FontSize',9,'BackgroundColor','w','Color','k','Clipping','on');
            end
        end
    end
    function patch_boundaries(S)
        t = linspace(0,1,20); uv = [t 1-t 0*t;0*t t 1-t];
        pol = surfsmooth3d.internal.koorn.pols(double(S.norders(1)),uv); xyz = nan(3,size(uv,2)+1,S.npatches);
        for k = 1:S.npatches, xyz(:,1:end-1,k) = S.srccoefs{k}(1:3,:)*pol; end
        xyz = reshape(xyz,3,[]); plot3(ax,xyz(1,:),xyz(2,:),xyz(3,:),'Color',[.2 .2 .2],'LineWidth',.4);
    end
    function light_update()
        delete(findobj(ax,'Type','light')); camlight(ax,'headlight'); lighting(ax,'gouraud');
    end
    function update_status()
        if ~isgraphics(fig), return; end
        message = state.notice;
        if ~isempty(state.result)
            r = state.result; info = r.info; health = 'PASS';
            if r.report.has_warnings, health = 'WARN'; end
            if ~r.report.is_valid, health = 'FAIL'; end
            message = sprintf('%s\nDrawn: p=%d, rlam=%.4g, %s, smooth edges %s\n%d patches, %d nodes; depth %d; unresolved %d; partial validity %s', ...
                message,r.settings.order,r.settings.rlam,surfsmooth3d.edgepreserve.sigma_mode_label(r.settings.adapt_sigma), ...
                mat2str(r.settings.selectedEdgeIds),r.surface.npatches,r.surface.npts,info.achieved_depth,info.unresolved_count,health);
            message = sprintf('%s\nPartial indicators [position tail, normal tail, position check, normal check]: %s',message,mat2str(max(info.indicators,[],1),3));
            try
                different = ~isequal(read_settings(),r.settings);
            catch
                different = true;
            end
            if different, message = [message newline 'Pending settings: recompute to assess them. Save uses the previous drawn settings.']; end
        end
        status.String = message;
        controls.draw.Enable = onoff(~state.busy && ~isempty(state.pkg));
        controls.save.Enable = onoff(~state.busy && ~isempty(state.result));
    end
    function close_viewer(varargin)
        if state.busy, return; end
        delete(fig);
    end
end

function value = onoff(flag)
if flag, value = 'on'; else, value = 'off'; end
end

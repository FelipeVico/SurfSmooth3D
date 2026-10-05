function fig = launch_gui(repoRoot, packageDir, visibility, outputParent)
%LAUNCH_GUI Fixed-source preview and batch export using ordinary MATLAB UI.
if nargin < 2, packageDir = ''; end
if nargin < 3, visibility = 'on'; end
settings = surfsmooth3d.edgepreserve.default_settings();
sigmaFlags = [0 3 1];
state = struct('pkg',[],'work',[],'cleanup',[],'result',[], ...
    'drawn',[],'busy',false,'limits',[],'notice','Choose a geometry package.', ...
    'failureFigure',[]);
inputRoot = surfsmooth3d.edgepreserve.default_input_root(repoRoot);
if nargin < 4
    outputParent = fullfile(repoRoot,'outputs','smoother');
end
fig = figure('Name','Fixed-CAD two-stage smoothing','NumberTitle','off', ...
    'Position',[100 70 1440 960],'Visible',visibility, ...
    'CloseRequestFcn',@close_viewer);
ax = axes('Parent',fig,'Units','normalized','Position',[.05 .16 .53 .76]);
cb = colorbar(ax);
set(cb,'Units','normalized','Position',[.605 .16 .017 .76]);
cb.Ruler.TickLabelFormat = '%.3g';
panel = uipanel(fig,'Title','Fixed CAD source / two-stage smoother', ...
    'Units','normalized','Position',[.665 .025 .325 .95]);
status = uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[.035 .015 .595 .115],'HorizontalAlignment','left','FontSize',11,'Tag','status');
controls = struct();
button('Browse Package', [.04 .945 .44 .037],@browse_package,'browse');
button('Output Folder', [.52 .945 .44 .037],@choose_output,'output');
controls.source = uicontrol(panel,'Style','text','Units','normalized', ...
    'Position',[.04 .865 .92 .073],'HorizontalAlignment','left', ...
    'String','No package loaded.');
label('Smooth edge IDs', [.04 .835 .92 .025]);
controls.edges = edit_box('[]',[.04 .801 .92 .035],@cheap_changed,'edges');
controls.labels = checkbox('Show CAD edge labels',true, ...
    [.04 .768 .92 .028],@view_changed,'labels');
parameter('rlam',settings.rlam,[.5 20],.708,false);
label('Sigma mode', [.04 .668 .20 .03]);
controls.adapt = uicontrol(panel,'Style','popupmenu','Units','normalized', ...
    'Position',[.26 .668 .70 .032], ...
    'String',arrayfun(@surfsmooth3d.edgepreserve.sigma_mode_label,sigmaFlags,'UniformOutput',false), ...
    'Value',find(sigmaFlags==settings.adapt_sigma,1), ...
    'Callback',@expensive_changed,'Tag','adapt');
parameter('c_eps',settings.c_eps,[.05 3],.613,true);
parameter('c_rad',settings.c_rad,[.1 10],.558,true);
parameter('c_tau',settings.c_tau,[.05 5],.503,true);
label('Preview order / refinement', [.04 .472 .92 .025]);
controls.previewOrder = edit_box('6',[.04 .439 .43 .033],@expensive_changed,'previewOrder');
controls.previewLevel = edit_box('0',[.53 .439 .43 .033],@expensive_changed,'previewLevel');
label('Batch orders', [.04 .413 .92 .024]);
controls.orders = edit_box('[4 6 8]',[.04 .382 .92 .032],@expensive_changed,'orders');
label('Batch refinement levels', [.04 .357 .92 .024]);
controls.levels = edit_box('[1 2 3]',[.04 .326 .92 .032],@expensive_changed,'levels');
button('Draw / Recompute',[.04 .276 .44 .039],@draw_preview,'draw');
button('Save Current Pair',[.52 .276 .44 .039],@save_preview,'save');
button('Generate And Save Batch',[.04 .229 .92 .039],@generate_batch,'batch');
label('Display mode', [.04 .199 .92 .025]);
controls.display = uicontrol(panel,'Style','popupmenu','Units','normalized', ...
    'Position',[.04 .169 .92 .032],'String', ...
    {'partially smoothed','fully smooth','flat scaffold','CAD source'}, ...
    'Callback',@view_changed,'Tag','display');
label('Color mode', [.04 .140 .92 .025]);
controls.color = uicontrol(panel,'Style','popupmenu','Units','normalized', ...
    'Position',[.04 .110 .92 .032],'String', ...
    {'beta','signed mean curvature','log10 abs mean curvature', ...
     'patch max abs curvature','jacobian ratio to CAD','displacement to CAD','constant'}, ...
    'Value',3,'Callback',@view_changed,'Tag','color');
controls.wire = checkbox('Show flat scaffold wireframe',false, ...
    [.04 .074 .92 .027],@view_changed,'wire');
controls.showEdges = checkbox('Show CAD edges',true, ...
    [.04 .041 .92 .027],@view_changed,'showEdges');
rotate = rotate3d(fig);
rotate.Enable = 'on';
rotate.ActionPostCallback = @(~,~) update_light();
if ~isempty(packageDir)
    load_source(packageDir);
elseif strcmp(visibility,'on')
    browse_package();
else
    update_status();
end

    function h = label(text,position)
        h = uicontrol(panel,'Style','text','String',text,'Units','normalized', ...
            'Position',position,'HorizontalAlignment','left');
    end
    function h = button(text,position,callback,tag)
        h = uicontrol(panel,'Style','pushbutton','String',text,'Units','normalized', ...
            'Position',position,'Callback',callback,'Tag',tag);
    end
    function h = edit_box(text,position,callback,tag)
        h = uicontrol(panel,'Style','edit','String',text,'Units','normalized', ...
            'Position',position,'Callback',callback,'Tag',tag,'HorizontalAlignment','left');
    end
    function h = checkbox(text,value,position,callback,tag)
        h = uicontrol(panel,'Style','checkbox','String',text,'Value',value, ...
            'Units','normalized','Position',position,'Callback',callback,'Tag',tag);
    end
    function parameter(name,value,range,y,cheap)
        label(name,[.04 y+.031 .92 .022]);
        controls.(name) = edit_box(num2str(value),[.71 y .25 .032], ...
            @(~,~) parameter_edit(name,cheap),name);
        controls.([name 'Slider']) = uicontrol(panel,'Style','slider', ...
            'Units','normalized','Position',[.04 y .64 .03], ...
            'Min',range(1),'Max',range(2),'Value',value,'Tag',[name 'Slider'], ...
            'Callback',@(h,~) parameter_slide(h,name,cheap));
    end
    function parameter_edit(name,cheap)
        try
            value = positive_value(controls.(name));
            slider = controls.([name 'Slider']);
            set(slider,'Min',min(slider.Min,value),'Max',max(slider.Max,value),'Value',value);
            if cheap, cheap_changed(); else, expensive_changed(); end
        catch exception
            show_error(exception);
        end
    end
    function parameter_slide(h,name,cheap)
        set(controls.(name),'String',sprintf('%.6g',h.Value));
        if cheap, cheap_changed(); else, expensive_changed(); end
    end
    function current = read_controls()
        current = settings;
        current.selectedEdgeIds = surfsmooth3d.edgepreserve.parse_integer_list(controls.edges.String,1,inf,true);
        if ~isempty(state.pkg) && any(~ismember(current.selectedEdgeIds,[state.pkg.edges.edge_id]))
            error('edgepreserve:edgeSelection','Selected CAD edge IDs are not in this package.');
        end
        for name = {'rlam','c_eps','c_rad','c_tau'}
            current.(name{1}) = positive_value(controls.(name{1}));
        end
        current.adapt_sigma = sigmaFlags(controls.adapt.Value);
        current.preview_order = single_integer(controls.previewOrder.String,1,20);
        current.preview_level = single_integer(controls.previewLevel.String,0,inf);
        current.orders = surfsmooth3d.edgepreserve.parse_integer_list(controls.orders.String,1,20);
        current.levels = surfsmooth3d.edgepreserve.parse_integer_list(controls.levels.String,0,inf);
    end
    function browse_package(varargin)
        start = inputRoot;
        if ~isempty(state.pkg), start = state.pkg.package_dir; end
        if ~isfolder(start), start = fileparts(repoRoot); end
        folder = uigetdir(start,'Select a geometry package folder');
        if isequal(folder,0), return; end
        load_source(folder);
    end
    function load_source(folder)
        try
            pkg = surfsmooth3d.edgepreserve.load_package(folder);
            state.cleanup = [];
            state.work = [];
            state.pkg = pkg;
            state.result = [];
            state.drawn = [];
            controls.edges.String = '[]';
            controls.source.String = sprintf('%s / %s\nCAD order %d, source refinement %d\nCoarse triangles %d, CAD edges %d', ...
                pkg.geometry_name,pkg.name,pkg.source_order,pkg.source_level, ...
                pkg.coarse_triangle_count,numel(pkg.edges));
            state.limits = limits_for({pkg.surface.r,pkg.scaffold.nodes});
            state.notice = 'CAD loaded. Select edges to smooth, then Draw / Recompute.';
            cla(ax);
            refresh_plot();
        catch exception
            show_error(exception);
        end
    end
    function choose_output(varargin)
        start = outputParent;
        if ~isfolder(start), start = repoRoot; end
        folder = uigetdir(start,'Choose the output parent folder');
        if ~isequal(folder,0)
            outputParent = folder;
            state.notice = ['Output parent: ' folder];
            update_status();
        end
    end
    function expensive_changed(varargin)
        if state.busy, return; end
        try
            read_controls();
            state.notice = 'Settings edited. The smoother runs only on Draw / Recompute or batch launch.';
            update_status();
        catch exception
            show_error(exception);
        end
    end
    function cheap_changed(varargin)
        if state.busy, return; end
        try
            current = read_controls();
            if ~isempty(state.result)
                state.result = surfsmooth3d.edgepreserve.update_blend(state.result,state.pkg.edges,current);
            end
            state.notice = 'Blend controls updated; fully smooth surface unchanged.';
            refresh_plot();
        catch exception
            show_error(exception);
        end
    end
    function view_changed(varargin)
        if ~state.busy, refresh_plot(); end
    end
    function draw_preview(varargin)
        if state.busy || isempty(state.pkg), return; end
        try
            current = read_controls();
            if ~confirm_size(current.preview_order,current.preview_level), return; end
            set_busy(true);
            release = onCleanup(@() set_busy(false));
            progress('Running two-stage preview...');
            if isempty(state.work)
                [state.work,state.cleanup] = surfsmooth3d.edgepreserve.stage_source(state.pkg);
            end
            opts = surfsmooth3d.edgepreserve.mesher_options(state.pkg,state.work,current,current.preview_level);
            all = surfsmooth3d.multiscale_mesher(state.work.scaffold_file,current.preview_order,opts);
            result = surfsmooth3d.edgepreserve.compute_pair(state.pkg,state.work,all{end}, ...
                current.preview_level,current);
            state.result = result;
            state.drawn = current;
            state.limits = limits_for({state.pkg.surface.r,state.pkg.scaffold.nodes, ...
                result.S_smooth.r,result.S_cad.r});
            surfsmooth3d.edgepreserve.print_validity_report(result.report_smooth);
            if ~isempty(result.report_blend), surfsmooth3d.edgepreserve.print_validity_report(result.report_blend); end
            state.notice = 'Preview drawn. Both surfaces can be saved independently if valid.';
            refresh_plot();
        catch exception
            show_error(exception);
        end
    end
    function save_preview(varargin)
        if state.busy || isempty(state.result), return; end
        try
            % Save exactly the latest drawn expensive settings and live blend.
            current = read_controls();
            state.result = surfsmooth3d.edgepreserve.update_blend(state.result,state.pkg.edges,current);
            snapshot = state.result.params;
            snapshot = rmfield(snapshot,{'sigma','geometryDiameter'});
            snapshot.orders = state.result.order;
            snapshot.levels = state.result.level;
            set_busy(true);
            release = onCleanup(@() set_busy(false));
            directory = surfsmooth3d.edgepreserve.new_run_directory(outputParent,state.pkg,snapshot,'preview');
            rows = surfsmooth3d.edgepreserve.save_pair(directory,state.pkg,state.result);
            manifest = struct2table(rows);
            writetable(manifest,fullfile(directory,'batch_manifest.csv'));
            state.notice = sprintf('Saved %d/2 displayed-preview surfaces: %s', ...
                nnz(manifest.status == "saved"),directory);
            fprintf('%s\n',state.notice);
            disp(manifest(:,{'mode','status','error_message'}));
            refresh_plot();
        catch exception
            show_error(exception);
        end
    end
    function generate_batch(varargin)
        if state.busy || isempty(state.pkg), return; end
        try
            snapshot = read_controls();
            if ~confirm_size(max(snapshot.orders),max(snapshot.levels)), return; end
            set_busy(true);
            release = onCleanup(@() set_busy(false));
            [manifest,directory] = surfsmooth3d.edgepreserve.run_batch( ...
                state.pkg,snapshot,outputParent,@progress,struct('failure',@show_newton_failure));
            state.notice = sprintf('Batch saved %d/%d surfaces: %s', ...
                nnz(manifest.status == "saved"),height(manifest),directory);
            update_status();
        catch exception
            show_error(exception);
        end
    end
    function refresh_plot()
        if isempty(state.pkg), update_status(); return; end
        camera = struct('position',ax.CameraPosition,'target',ax.CameraTarget, ...
            'up',ax.CameraUpVector,'angle',ax.CameraViewAngle);
        hadPlot = ~isempty(ax.Children);
        cla(ax);
        hold(ax,'on');
        mode = controls.display.Value;
        colorNames = controls.color.String;
        colorMode = colorNames{controls.color.Value};
        beta = [];
        report = [];
        note = '';
        if mode == 3
            surfsmooth3d.edgepreserve.plot_gidmsh(ax,state.pkg.scaffold,'solid');
            values = zeros(1,1); colorMode = 'constant'; colorLabel = 'flat scaffold';
            titleText = 'flat scaffold: no high-order curvature';
        else
            S = state.pkg.surface;
            reference = S;
            titleText = 'fixed CAD source (preview not computed)';
            if mode == 4
                titleText = 'fixed CAD source';
            elseif ~isempty(state.result)
                reference = state.result.S_cad;
                if mode == 1 && ~isempty(state.result.S_blend)
                    S = state.result.S_blend; beta = state.result.beta;
                    report = state.result.report_blend;
                    titleText = 'partially smoothed';
                else
                    S = state.result.S_smooth; report = state.result.report_smooth;
                    titleText = 'fully smooth';
                    if mode == 1, note = ['Partial build failed: ' state.result.blend_error]; end
                end
            end
            [values,colorLabel,colorNote] = surfsmooth3d.edgepreserve.color_values(S,reference,beta,colorMode);
            note = strjoin({note,colorNote},' ');
            axes(ax);
            surfsmooth3d.edgepreserve.plot_surface(S,values,{},max(8,2*double(max(S.norders))));
            if controls.wire.Value, surfsmooth3d.edgepreserve.plot_gidmsh(ax,state.pkg.scaffold,'wireframe'); end
        end
        surfsmooth3d.edgepreserve.color_limits(ax,values,colorMode,strtrim(note));
        cb.Label.String = colorLabel;
        set(cb,'Position',[.605 .16 .017 .76]);
        set(ax,'Position',[.05 .16 .53 .76]);
        plot_edges();
        axis(ax,'equal');
        xlim(ax,state.limits(1,:)); ylim(ax,state.limits(2,:)); zlim(ax,state.limits(3,:));
        axis(ax,'vis3d'); grid(ax,'on');
        xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
        if hadPlot
            set(ax,'CameraPosition',camera.position,'CameraTarget',camera.target, ...
                'CameraUpVector',camera.up,'CameraViewAngle',camera.angle);
        else
            view(ax,3); ax.CameraTarget = mean(state.limits,2).';
        end
        title(ax,{sprintf('%s: %s',state.pkg.geometry_name,titleText),colorLabel}, ...
            'Interpreter','none','FontSize',12);
        update_light();
        lighting(ax,'gouraud');
        setappdata(fig,'edgepreserve_preview',state.result);
        setappdata(fig,'edgepreserve_source',state.pkg);
        update_status(report,strtrim(note));
        drawnow limitrate;
    end
    function plot_edges()
        if ~controls.showEdges.Value, return; end
        selected = [];
        try
            current = read_controls();
            selected = current.selectedEdgeIds;
        catch
            % Invalid pending IDs must not prevent viewing the drawn geometry.
        end
        for edge = state.pkg.edges(:).'
            color = [.4 .4 .4]; width = .7;
            if ismember(edge.edge_id,selected), color = [0 0 0]; width = 2; end
            xyz = [];
            for piece = edge.panels(:).'
                line(ax,piece.xyz(1,:),piece.xyz(2,:),piece.xyz(3,:), ...
                    'Color',color,'LineWidth',width,'HitTest','off');
                xyz = [xyz piece.xyz]; %#ok<AGROW>
            end
            if controls.labels.Value && ~isempty(xyz)
                point = edge_midpoint(edge);
                text(ax,point(1),point(2),point(3),sprintf('E%d',edge.edge_id), ...
                    'BackgroundColor','w','Color','k','FontSize',9,'HitTest','off');
            end
        end
    end
    function update_light()
        delete(findall(ax,'Type','light'));
        camlight(ax,'headlight');
    end
    function update_status(report,note)
        if nargin < 1, report = []; end
        if nargin < 2, note = ''; end
        pending = '';
        if ~isempty(state.result)
            try
                current = read_controls();
                names = {'rlam','adapt_sigma','preview_order','preview_level'};
                mismatch = names(cellfun(@(key) ~isequal(current.(key),state.drawn.(key)),names));
                if isempty(mismatch)
                    pending = 'Drawn preview matches expensive controls.';
                else
                    pending = ['PENDING (not used in current pair): ' strjoin(mismatch,', ')];
                end
            catch exception
                pending = exception.message;
            end
            drawn = sprintf('Drawn: p%d, refinement %d, rlam %.4g, %s, two-stage ON', ...
                state.result.order,state.result.level,state.result.params.rlam, ...
                surfsmooth3d.edgepreserve.sigma_mode_label(state.result.params.adapt_sigma));
            sigma = state.result.sigma;
            if isempty(sigma)
                drawn = [drawn ' | sigma unavailable'];
            else
                drawn = sprintf('%s | sigma min/mean/max %.3g %.3g %.3g',drawn,min(sigma),mean(sigma),max(sigma));
            end
        else
            drawn = 'No smoother computed yet. Empty edge selection means CAD-only partial output.';
        end
        health = '';
        if ~isempty(report)
            health = health_text(report);
        end
        if ~isempty(state.result) && ~isempty(state.result.beta)
            beta = state.result.beta;
            health = sprintf('%s\nc_eps %.3g, c_rad %.3g, c_tau %.3g | beta min/mean/max %.3g %.3g %.3g', ...
                health,state.result.params.c_eps,state.result.params.c_rad, ...
                state.result.params.c_tau,min(beta),mean(beta),max(beta));
        end
        status.String = strjoin({state.notice,drawn,pending,health,note},newline);
        update_actions();
    end
    function progress(message)
        state.notice = message;
        fprintf('%s\n',message);
        update_status(); drawnow;
    end
    function show_error(exception)
        state.notice = ['Error: ' exception.message];
        if ismember(exception.identifier,{'MULTISCALE_MESHER:NewtonRadiusGuard', ...
                'MULTISCALE_MESHER:NewtonProjectionFailure'})
            report = surfsmooth3d.edgepreserve.capture_newton_failure(exception);
            if ~isempty(report)
                show_newton_failure(report,report.filename);
            else
                state.notice = [state.notice newline 'Localization data unavailable; the solver stopped safely.'];
            end
            if ~isempty(state.result)
                state.notice = [state.notice newline ...
                    'Previous successful preview retained; failed settings were not applied.'];
            end
        end
        warning('edgepreserve:gui','%s',exception.message);
        update_status();
    end
    function show_newton_failure(report,filename)
        try
            state.failureFigure = surfsmooth3d.edgepreserve.plot_newton_failure( ...
                state.pkg,report,state.failureFigure,visibility);
            setappdata(fig,'edgepreserve_newton_failure',report);
            fprintf('Newton diagnostic: %s\n',filename);
        catch diagnosticError
            warning('edgepreserve:newtonDiagnostic','Cannot plot Newton failure: %s',diagnosticError.message);
        end
    end
    function set_busy(value)
        state.busy = value;
        if ~isgraphics(fig), return; end
        handles = findall(panel,'Type','uicontrol');
        if value, set(handles,'Enable','off'); set(fig,'Pointer','watch');
        else, set(handles,'Enable','on'); set(fig,'Pointer','arrow'); end
        update_actions();
        drawnow;
    end
    function update_actions()
        if state.busy, return; end
        sourceEnable = 'on';
        if isempty(state.pkg)
            sourceEnable = 'off';
        end
        set(findobj(panel,'Tag','draw'),'Enable',sourceEnable);
        set(findobj(panel,'Tag','batch'),'Enable',sourceEnable);
        saveEnable = 'on';
        if isempty(state.result)
            saveEnable = 'off';
        end
        set(findobj(panel,'Tag','save'),'Enable',saveEnable);
    end
    function ok = confirm_size(order,level)
        nodes = state.pkg.coarse_triangle_count*4^level*(order+1)*(order+2)/2;
        ok = true;
        if nodes > 2e6
            answer = questdlg(sprintf('Largest surface has %.3g nodes. Continue?',nodes), ...
                'Large resolution study','Continue','Cancel','Cancel');
            ok = strcmp(answer,'Continue');
        end
    end
    function close_viewer(varargin)
        if state.busy
            state.notice = 'A computation is active. Close the viewer after it completes.';
            update_status(); return;
        end
        state.cleanup = [];
        delete(fig);
    end
end

function value = positive_value(control)
value = str2double(control.String);
if ~isscalar(value) || ~isfinite(value) || value<=0
    error('edgepreserve:positiveParameter', '%s must be positive and finite.',control.Tag);
end
end

function value = single_integer(text,lo,hi)
value = surfsmooth3d.edgepreserve.parse_integer_list(text,lo,hi);
if numel(value) ~= 1, error('edgepreserve:singleInteger','Preview resolution needs one integer.'); end
end

function limits = limits_for(arrays)
lo = inf(3,1); hi = -inf(3,1);
for k = 1:numel(arrays)
    lo = min(lo,min(arrays{k},[],2)); hi = max(hi,max(arrays{k},[],2));
end
center = (lo+hi)/2;
width = max(max(hi-lo),1)*1.12;
limits = [center-width/2 center+width/2];
end

function point = edge_midpoint(edge)
lengths = arrayfun(@(panel) sum(panel.w),edge.panels);
target = sum(lengths)/2;
index = find(cumsum(lengths)>=target,1);
if isempty(index), index = 1; end
previous = sum(lengths(1:index-1));
panel = edge.panels(index);
[~,node] = min(abs(cumsum(panel.w)-target+previous));
point = panel.xyz(:,node);
end

function text = health_text(report)
status = 'PASS';
if ~report.is_valid, status = 'FAIL'; elseif report.has_warnings, status = 'WARN'; end
values = accumarray(report.patch_id,abs(report.mean_curv),[],@max);
[maximum,patch] = max(values);
text = sprintf('Validity %s | max |H| %.3g (patch %d) | min jac/ref %.3g', ...
    status,maximum,patch,report.stats.jac_ratio_min);
end

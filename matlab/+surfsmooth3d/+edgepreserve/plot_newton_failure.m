function fig = plot_newton_failure(pkg, report, fig, visibility)
%PLOT_NEWTON_FAILURE Color launch points, never the runaway target coordinates.
if nargin < 3, fig = []; end
if nargin < 4, visibility = 'on'; end
if isempty(fig) || ~isgraphics(fig,'figure')
    fig = figure('Name','Two-stage Newton failure','NumberTitle','off', ...
        'Position',[130 100 1100 800],'Visible',visibility);
else
    clf(fig);
end
ax = axes('Parent',fig,'Position',[.075 .16 .73 .74]);
axes(ax);
surfsmooth3d.edgepreserve.plot_surface(pkg.surface,zeros(1,pkg.surface.npts), ...
    {'FaceColor',[.65 .68 .70],'FaceAlpha',.22,'HandleVisibility','off'}, ...
    max(8,2*double(max(pkg.surface.norders))));
hold(ax,'on');
finiteBase = all(isfinite(report.base),1).';
active = ~report.converged;
colored = finiteBase & active & ~report.rejected & report.step_available & ...
    isfinite(report.step_length) & report.step_length >= 0;
converged = finiteBase & report.converged & ~report.rejected;
unavailable = finiteBase & active & ~report.rejected & ~colored;
rejected = finiteBase & report.rejected;
handles = gobjects(0);
if any(converged)
    handles(end+1) = points(converged,10,[.48 .48 .48],'.','Converged');
end
if any(colored)
    % Compute in the log domain so normalization cannot overflow or underflow.
    colors = max(log10(report.step_length(colored))-log10(report.radius),-16);
    handles(end+1) = points(colored,22,colors,'o','Active: latest Newton step');
    limits = [min(colors) max(colors)];
    if limits(1)==limits(2), limits = limits + [-.5 .5]; end
else
    limits = [-16 -15];
end
if any(unavailable)
    handles(end+1) = points(unavailable,28,[.2 .2 .2],'d','Step unavailable / nonfinite');
end
if any(rejected)
    label = 'Rejected: outside 10R / nonfinite';
    if strcmp(report.failure_kind,'newton'), label = 'Failed: Newton stopping criterion'; end
    handles(end+1) = points(rejected,70,[1 0 0],'o',label);
end
colormap(ax,parula);
clim(ax,limits);
cb = colorbar(ax);
cb.Label.String = 'log10(max(latest physical Newton step / R, 1e-16))';
if ~any(colored), cb.Label.String = [cb.Label.String ' (no finite active step data)']; end
if ~isempty(handles), legend(ax,handles,'Location','best'); end
original = [pkg.surface.r pkg.scaffold.nodes];
lower = min(original,[],2); upper = max(original,[],2);
padding = .04*max(upper-lower);
if padding <= 0, padding = .04*report.radius; end
axis(ax,'equal');
xlim(ax,[lower(1)-padding upper(1)+padding]);
ylim(ax,[lower(2)-padding upper(2)+padding]);
zlim(ax,[lower(3)-padding upper(3)+padding]);
axis(ax,'vis3d'); view(ax,3); grid(ax,'on');
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
ax.CameraTarget = ((lower+upper)/2).';
title(ax,{sprintf('%s: Newton stopped',pkg.geometry_name), ...
    sprintf('%s, refinement %d, iteration %d; red = rejected launch points', ...
    report.stage,report.refinement,report.iteration),report.failure_reason},'Interpreter','none');
camlight(ax,'headlight'); lighting(ax,'gouraud'); rotate3d(fig,'on');
finiteSteps = report.step_length(report.step_available & isfinite(report.step_length));
stats = 'No finite step measurements available.';
if ~isempty(finiteSteps)
    stats = sprintf('Finite physical step min/median/max: %.4g / %.4g / %.4g', ...
        min(finiteSteps),median(finiteSteps),max(finiteSteps));
end
summary = sprintf(['Failed/rejected %d/%d; active %d; fixed R = %.4g, radius guard = %.4g\n' ...
    '%s\nLarge steps indicate suspicious regions, not confirmed failure.'], ...
    nnz(report.rejected),numel(report.index),nnz(active),report.radius, ...
    report.limit_multiplier*report.radius,stats);
uicontrol(fig,'Style','text','Units','normalized','Position',[.075 .02 .86 .10], ...
    'String',summary,'HorizontalAlignment','left','FontSize',11);
fprintf('Newton failure: %s, refinement %d, iteration %d\n%s\n', ...
    report.stage,report.refinement,report.iteration,summary);
indices = find(report.rejected);
fprintf('Rejected point indices and distance/R (first %d):\n',min(30,numel(indices)));
fprintf('  %d  %.6g\n',[report.index(indices(1:min(30,end))) ...
    report.distance_over_radius(indices(1:min(30,end)))].');
setappdata(fig,'edgepreserve_newton_failure',report);
drawnow;

    function h = points(mask,size,color,marker,name)
        xyz = report.base(:,mask);
        h = scatter3(ax,xyz(1,:),xyz(2,:),xyz(3,:),size,color,marker,'filled', ...
            'DisplayName',name,'Tag',name);
    end
end

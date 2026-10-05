% Fixed-CAD two-stage smoother: interactive preview and resolution study.
% The selected package supplies ONE fixed CAD quadrature, original scaffold,
% and edge quadrature. Output order/refinement are independent GUI choices.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot,'matlab')); setup_surfsmooth3d;

% Leave empty to choose a package folder when the viewer opens.
packageDir = '';
% packageRoot = surfsmooth3d.edgepreserve.default_input_root(repoRoot);
% packageDir = fullfile(packageRoot, ...
%     'convergence_test_geometry_manas','p8_mf0p10000_rl0');
% packageDir = fullfile(packageRoot, ...
%     'convergence_intersection_two_balls','p8_mf0p10000_rl1');

viewer = surfsmooth3d.edgepreserve.launch_gui(repoRoot,packageDir);

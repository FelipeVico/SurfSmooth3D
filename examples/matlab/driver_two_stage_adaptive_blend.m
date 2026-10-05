% Adaptively refine the partially smoothed surface; selected edges are smoothed.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot,'matlab')); setup_surfsmooth3d;
packageDir = ''; % Select the fixed CAD geometry package in the GUI.
maxPoints = 2000000; % Advanced cap, checked before each child solve batch.
viewer = surfsmooth3d.adaptiveblend.launch_gui(repoRoot,packageDir,'on',[],maxPoints);

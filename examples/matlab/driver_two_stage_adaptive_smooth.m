% Fully smooth, fixed-order adaptive refinement of one fixed CAD level set.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot,'matlab')); setup_surfsmooth3d;
packageDir = ''; % Choose a geometry package folder in the GUI.
% Advanced resource cap, applied before allocating a child solve batch.
maxPoints = 2000000;
viewer = surfsmooth3d.adaptivesmoother.launch_gui(repoRoot,packageDir,'on',[],maxPoints);

function root = setup_surfsmooth3d()
%SETUP_SURFSMOOTH3D Add only this installation's MATLAB package root.
matlabRoot = fileparts(mfilename('fullpath'));
addpath(matlabRoot);
root = fileparts(matlabRoot);
end

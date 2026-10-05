function h = plot_surface(mesh, varargin)
%PLOT_SURFACE Plot the surfer surface using SurfSmooth3D geometry utilities.

S = surfsmooth3d.stepmesher.to_surfer(mesh);
h = plot(S, varargin{:});
end

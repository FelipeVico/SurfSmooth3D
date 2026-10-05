function h = plot_projection_error(mesh, varargin)
%PLOT_PROJECTION_ERROR Scatter projected nodes colored by projection distance.

xyz = double(mesh.xyz);
dist = double(mesh.proj_dist);
h = scatter3(xyz(1, :), xyz(2, :), xyz(3, :), 16, dist, 'filled', varargin{:});
axis equal
colorbar
xlabel('x'); ylabel('y'); zlabel('z');
end

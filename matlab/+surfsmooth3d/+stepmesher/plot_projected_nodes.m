function h = plot_projected_nodes(mesh, varargin)
%PLOT_PROJECTED_NODES Scatter plot of projected RV nodes.

xyz = double(mesh.xyz);
h = scatter3(xyz(1, :), xyz(2, :), xyz(3, :), 12, varargin{:});
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
end

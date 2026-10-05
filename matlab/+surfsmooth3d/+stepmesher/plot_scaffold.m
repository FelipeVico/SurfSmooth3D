function h = plot_scaffold(mesh, varargin)
%PLOT_SCAFFOLD Plot the linear Gmsh surface scaffold.

nodes = double(mesh.scaffold_nodes);
tris = double(mesh.scaffold_triangles).';
if isempty(nodes) || isempty(tris)
    error('stepmesher:badMesh', 'Mesh does not contain scaffold data.');
end
h = trisurf(tris, nodes(1, :), nodes(2, :), nodes(3, :), varargin{:});
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
end

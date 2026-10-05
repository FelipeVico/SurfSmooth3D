function mesh = mesh_step(stepFile, opts)
%MESH_STEP Mesh a STEP file and return projected RV triangle nodes.
%
% mesh = surfsmooth3d.stepmesher.mesh_step(stepFile, opts)
%
% The MEX gateway uses the C++ core for STEP import, Gmsh/OCC scaffolding,
% and CAD surface evaluation. MATLAB supplies the exact RV nodes from
% fmm3dbie's koorn package so the returned samples are ready for surfsmooth3d.surfer.

if nargin < 2 || isempty(opts)
    opts = struct();
end
opts = default_options(opts);

if isempty(which('surfsmooth3d.internal.koorn.rv_nodes'))
    error('stepmesher:missingKoorn', ...
        ['surfsmooth3d.internal.koorn.rv_nodes was not found. Add SurfSmooth3D/matlab to the ', ...
         'MATLAB path before calling surfsmooth3d.stepmesher.mesh_step.']);
end

uv = surfsmooth3d.internal.koorn.rv_nodes(opts.order);
if isempty(which('step_mesher_mex'))
    error('stepmesher:missingMex', ...
        ['step_mesher_mex is not built. Run CMake with ', ...
         'STEP_MESHER_BUILD_MEX=ON, then add this repo''s matlab folder ', ...
         'to the MATLAB path.']);
end

mesh = step_mesher_mex('mesh_step', char(stepFile), opts, uv);
mesh.uv = uv;
mesh.nodes_per_patch = double(mesh.nodes_per_patch);
mesh.npatches = double(mesh.npatches);
mesh.norder = double(mesh.norder);
mesh.iptype = double(mesh.iptype);
mesh.edge_order = double(mesh.edge_order);
mesh.edge_nodes_per_panel = double(mesh.edge_nodes_per_panel);
mesh.ncad_edges = double(mesh.ncad_edges);
mesh.npanels_total = double(mesh.npanels_total);
end

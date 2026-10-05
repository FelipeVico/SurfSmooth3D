function smoke_mesh_step(repoRoot)
%SMOKE_MESH_STEP Basic MATLAB smoke test.

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
addpath(fullfile(repoRoot, 'matlab'));

stepFile = fullfile(repoRoot, 'step_mesher', 'examples', 'step_files', 'test_geometry_1.step');
opts = struct('order', 2, 'mesh_fraction', 0.20, 'refinement_level', 0);
mesh = surfsmooth3d.stepmesher.mesh_step(stepFile, opts);

assert(isfield(mesh, 'xyz'));
assert(size(mesh.xyz, 1) == 3);
assert(mesh.npatches > 0);
assert_koorn_vertex_convention(mesh, opts.order);
assert_quality_fields(mesh);

optsCurvature = opts;
optsCurvature.mesh_profile = 'curvature';
optsCurvature.curvature_points = 20;
optsCurvature.hmin_fraction = 0.01;
optsCurvature.hmax_fraction = 0.25;
optsCurvature.gmsh_algorithm = 'frontal_delaunay';
optsCurvature.optimize = true;
meshCurvature = surfsmooth3d.stepmesher.mesh_step(stepFile, optsCurvature);
assert(strcmp(meshCurvature.meta.mesh_profile, 'curvature'));
assert(meshCurvature.meta.mesh_size_extend_from_boundary);
assert_quality_fields(meshCurvature);

[optsLocal, localChoice] = surfsmooth3d.stepmesher.apply_mesh_profile_choice( ...
    opts, 'curvature_balanced_local');
assert(strcmp(localChoice.name, "curvature_balanced_local"));
assert(~optsLocal.mesh_size_extend_from_boundary);
meshLocal = surfsmooth3d.stepmesher.mesh_step(stepFile, optsLocal);
assert(strcmp(meshLocal.meta.mesh_profile, 'curvature'));
assert(~meshLocal.meta.mesh_size_extend_from_boundary);
assert_quality_fields(meshLocal);

[optsCadEdge, cadEdgeChoice] = surfsmooth3d.stepmesher.apply_mesh_profile_choice( ...
    opts, 'curvature_balanced_local_cadedge');
assert(strcmp(cadEdgeChoice.name, "curvature_balanced_local_cadedge"));
assert(strcmp(optsCadEdge.hmin_mode, 'cad_edge'));
meshCadEdge = surfsmooth3d.stepmesher.mesh_step(stepFile, optsCadEdge);
assert(strcmp(meshCadEdge.meta.hmin_mode, 'cad_edge'));
assert(meshCadEdge.meta.meaningful_cad_edge_count > 0);
expectedHmin = meshCadEdge.meta.hmin_edge_scale * ...
    meshCadEdge.meta.min_meaningful_cad_edge_length;
assert(abs(meshCadEdge.meta.hmin_effective - expectedHmin) <= ...
    1e-12 * max(1, abs(expectedHmin)));
assert_quality_fields(meshCadEdge);

optsRefined = opts;
optsRefined.refinement_level = 1;
meshRefined = surfsmooth3d.stepmesher.mesh_step(stepFile, optsRefined);
assert_refined_koorn_vertex_convention(meshRefined, optsRefined.order, ...
    optsRefined.refinement_level);

[srcvals, norders, iptype] = surfsmooth3d.stepmesher.to_srcvals(mesh);
assert(size(srcvals, 1) == 12);
assert(all(isfinite(srcvals(:))));
assert(all(norders == opts.order));
assert(all(iptype == 1));

stats = surfsmooth3d.stepmesher.validate_area_flux(mesh);
assert(isfinite(stats.area) && stats.area > 0);

go3File = [tempname, '.go3'];
gidmshFile = [tempname, '.gidmsh'];
edgesFile = [tempname, '.txt'];
packageDir = [tempname, '_geom_package'];
cleanup = onCleanup(@() cleanup_paths({go3File, gidmshFile, edgesFile, ...
    packageDir}));
payload = surfsmooth3d.stepmesher.export_go3(mesh, go3File);
assert(isfile(go3File));
assert(payload.total_nodes == size(srcvals, 2));

S = surfsmooth3d.surfer.load_from_file(go3File);
tol = 1e-10;
r = srcvals(1:3, :);
du = srcvals(4:6, :);
dv = srcvals(7:9, :);
normal = srcvals(10:12, :);
assert(max(abs(S.r(:) - r(:))) < tol);
assert(max(abs(S.du(:) - du(:))) < tol);
assert(max(abs(S.dv(:) - dv(:))) < tol);
assert(max(abs(S.n(:) - normal(:))) < tol);

gidmshPayload = surfsmooth3d.stepmesher.export_gidmsh(mesh, gidmshFile);
assert(isfile(gidmshFile));
assert(gidmshPayload.num_nodes == size(mesh.scaffold_nodes, 2));
assert(gidmshPayload.num_triangles == size(mesh.scaffold_triangles, 2));
txt = fileread(gidmshFile);
assert(contains(txt, 'MESH dimension 3 ElemType Triangle Nnode 3'));
assert(contains(txt, 'Coordinates'));
assert(contains(txt, 'Elements'));

edgesPayload = surfsmooth3d.stepmesher.export_edges_gll(mesh, edgesFile);
assert(isfile(edgesFile));
assert(edgesPayload.num_cad_edges == mesh.ncad_edges);
assert(edgesPayload.num_panels == mesh.npanels_total);
edges = surfsmooth3d.stepmesher.load_edges_gll(edgesFile);
assert(numel(edges) == mesh.ncad_edges);
assert(all([edges.num_panels] >= 1));
assert(count_panels(edges) == mesh.npanels_total);
assert(all_panels_valid(edges, mesh.edge_nodes_per_panel));

packagePayload = surfsmooth3d.stepmesher.export_geom_package(mesh, packageDir, ...
    'smoke_mesh_step');
assert(isfolder(packageDir));
assert(isfile(packagePayload.surface_file));
assert(isfile(packagePayload.scaffold_file));
assert(isfile(packagePayload.edges_file));
assert(isfile(packagePayload.metadata_file));
pkg = surfsmooth3d.stepmesher.load_geom_package(packageDir);
assert(numel(pkg.edges) == mesh.ncad_edges);

oldFigureVisible = get(0, 'DefaultFigureVisible');
visibilityCleanup = onCleanup(@() set(0, 'DefaultFigureVisible', ...
    oldFigureVisible));
set(0, 'DefaultFigureVisible', 'off');
viewerOpts = struct('surfaceRefine', 2, 'edgeRefine', 8, ...
    'showEdgeLabels', false, 'showEdgeNodes', true);
[viewerHandles, viewerData] = surfsmooth3d.stepmesher.plot_geom_package(pkg, ...
    viewerOpts);
figureCleanup = onCleanup(@() close_figure(viewerHandles.figure));
assert(size(viewerData.surface_vertices, 2) == 3);
assert(size(viewerData.surface_faces, 2) == 3);
assert(all(isfinite(viewerData.surface_vertices(:))));
assert(numel(viewerHandles.edges) == mesh.ncad_edges);
assert(all(arrayfun(@(h) isgraphics(h), viewerHandles.edges)));
qualityHandles = surfsmooth3d.stepmesher.plot_scaffold_quality(mesh, 'min_angle');
assert(isgraphics(qualityHandles));
close(viewerHandles.figure);
clear figureCleanup
clear visibilityCleanup
end

function assert_quality_fields(mesh)
ntri = size(mesh.scaffold_triangles, 2);
assert(isfield(mesh, 'scaffold_min_angle_deg'));
assert(isfield(mesh, 'scaffold_edge_ratio'));
assert(isfield(mesh, 'scaffold_triangle_area'));
assert(numel(mesh.scaffold_min_angle_deg) == ntri);
assert(numel(mesh.scaffold_edge_ratio) == ntri);
assert(numel(mesh.scaffold_triangle_area) == ntri);
assert(all(isfinite(mesh.scaffold_min_angle_deg)));
assert(all(isfinite(mesh.scaffold_edge_ratio)));
assert(all(isfinite(mesh.scaffold_triangle_area)));
assert(all(mesh.scaffold_min_angle_deg >= 0));
assert(all(mesh.scaffold_edge_ratio >= 1));
assert(all(mesh.scaffold_triangle_area > 0));
assert(isfield(mesh.meta, 'min_scaffold_angle_deg'));
assert(isfield(mesh.meta, 'max_scaffold_edge_ratio'));
assert(isfield(mesh.meta, 'bad_min_angle_count'));
assert(isfield(mesh.meta, 'bad_edge_ratio_count'));
assert(isfield(mesh.meta, 'worst_quality_triangle'));
assert(isfield(mesh.meta, 'mesh_size_extend_from_boundary'));
assert(isfield(mesh.meta, 'hmin_mode'));
assert(isfield(mesh.meta, 'hmin_effective'));
assert(isfield(mesh.meta, 'hmax_effective'));
assert(isfield(mesh.meta, 'hmin_edge_scale'));
assert(isfield(mesh.meta, 'cad_edge_tol_fraction'));
assert(isfield(mesh.meta, 'min_raw_cad_edge_length'));
assert(isfield(mesh.meta, 'min_meaningful_cad_edge_length'));
assert(isfield(mesh.meta, 'cad_edge_count'));
assert(isfield(mesh.meta, 'meaningful_cad_edge_count'));
assert(mesh.meta.min_scaffold_angle_deg == min(mesh.scaffold_min_angle_deg));
assert(mesh.meta.max_scaffold_edge_ratio == max(mesh.scaffold_edge_ratio));
end

function cleanup_paths(paths)
for i = 1:numel(paths)
    if isfolder(paths{i})
        rmdir(paths{i}, 's');
    elseif isfile(paths{i})
        delete(paths{i});
    end
end
end

function close_figure(fig)
if isgraphics(fig)
    close(fig);
end
end

function n = count_panels(edges)
n = 0;
for i = 1:numel(edges)
    n = n + numel(edges(i).panels);
end
end

function ok = all_panels_valid(edges, nodesPerPanel)
ok = true;
for i = 1:numel(edges)
    for j = 1:numel(edges(i).panels)
        panel = edges(i).panels(j);
        ok = ok && panel.num_nodes == nodesPerPanel;
        ok = ok && all(isfinite(panel.xyz(:)));
        ok = ok && all(isfinite(panel.tangent(:)));
        ok = ok && all(isfinite(panel.w(:)));
        ok = ok && all(panel.w(:) >= 0);
        tangentNorm = vecnorm(panel.tangent, 2, 1);
        ok = ok && all(abs(tangentNorm - 1) < 1e-10);
    end
end
end

function assert_koorn_vertex_convention(mesh, order)
uv = surfsmooth3d.internal.koorn.rv_nodes(order);
A = surfsmooth3d.internal.koorn.vals2coefs(order, uv);
vertexUv = [0, 1, 0; 0, 0, 1];
vertexPols = surfsmooth3d.internal.koorn.pols(order, vertexUv);
nodesPerPatch = double(mesh.nodes_per_patch);
ncheck = min(double(mesh.npatches), 20);
errs = zeros(ncheck, 1);
for patch = 1:ncheck
    cols = (patch - 1) * nodesPerPatch + (1:nodesPerPatch);
    coef = double(mesh.xyz(:, cols)) * A.';
    patchVerts = coef * vertexPols;
    parentTri = double(mesh.parent_coarse_tri(patch));
    triNodes = double(mesh.scaffold_triangles(:, parentTri));
    expectedVerts = double(mesh.scaffold_nodes(:, triNodes));
    errs(patch) = max(vecnorm(patchVerts - expectedVerts, 2, 1));
end
scale = max(1.0, double(mesh.meta.bbox_diag));
tol = 1.0e-8 * scale;
assert(max(errs) < tol, ...
    ['RV node convention mismatch: extrapolated patch vertices do not ', ...
     'match scaffold vertices in fmm3dbie/koorn order. max err = %.6e'], ...
    max(errs));
end

function assert_refined_koorn_vertex_convention(mesh, order, refinementLevel)
uv = surfsmooth3d.internal.koorn.rv_nodes(order);
A = surfsmooth3d.internal.koorn.vals2coefs(order, uv);
vertexUv = [0, 1, 0; 0, 0, 1];
vertexPols = surfsmooth3d.internal.koorn.pols(order, vertexUv);
refTris = refined_reference_triangles(refinementLevel);
nref = size(refTris, 3);
nodesPerPatch = double(mesh.nodes_per_patch);
ncheck = min(double(mesh.npatches), max(20, 20 * nref));
errs = zeros(ncheck, 1);
for patch = 1:ncheck
    cols = (patch - 1) * nodesPerPatch + (1:nodesPerPatch);
    coef = double(mesh.xyz(:, cols)) * A.';
    patchVerts = coef * vertexPols;
    parentTri = double(mesh.parent_coarse_tri(patch));
    triNodes = double(mesh.scaffold_triangles(:, parentTri));
    parentVerts = double(mesh.scaffold_nodes(:, triNodes));
    refIndex = mod(patch - 1, nref) + 1;
    expectedVerts = parentVerts * refTris(:, :, refIndex);
    errs(patch) = max(vecnorm(patchVerts - expectedVerts, 2, 1));
end
scale = max(1.0, double(mesh.meta.bbox_diag));
tol = 1.0e-8 * scale;
assert(max(errs) < tol, ...
    ['Refined RV node convention mismatch: extrapolated subpatch ', ...
     'vertices do not match the refined scaffold vertices. max err = %.6e'], ...
    max(errs));
end

function refTris = refined_reference_triangles(refinementLevel)
n = 2 ^ refinementLevel;
refTris = zeros(3, 3, n * n);
itri = 0;
for i = 0:(n - 1)
    for j = 0:(n - i - 1)
        p00 = bary_from_uv(i / n, j / n);
        p10 = bary_from_uv((i + 1) / n, j / n);
        p01 = bary_from_uv(i / n, (j + 1) / n);
        itri = itri + 1;
        refTris(:, :, itri) = [p00, p10, p01];
        if i + j < n - 1
            p11 = bary_from_uv((i + 1) / n, (j + 1) / n);
            itri = itri + 1;
            refTris(:, :, itri) = [p10, p11, p01];
        end
    end
end
refTris = refTris(:, :, 1:itri);
end

function bary = bary_from_uv(u, v)
bary = [1 - u - v; u; v];
end

function stats = validate_area_flux(mesh)
%VALIDATE_AREA_FLUX Compute area and closed-surface normal flux diagnostics.

[srcvals] = surfsmooth3d.stepmesher.to_srcvals(mesh);
order = double(mesh.norder);
npatches = double(mesh.npatches);
nodes_per_patch = double(mesh.nodes_per_patch);
wref = surfsmooth3d.internal.koorn.rv_weights(order);
wref = wref(:).';

du = srcvals(4:6, :);
dv = srcvals(7:9, :);
normal_vec = cross(du, dv, 1);
jac = vecnorm(normal_vec, 2, 1);
normal = srcvals(10:12, :);

weights = zeros(1, size(srcvals, 2));
patchAreas = zeros(npatches, 1);
for patch = 1:npatches
    cols = (patch - 1) * nodes_per_patch + (1:nodes_per_patch);
    weights(cols) = jac(cols) .* wref;
    patchAreas(patch) = sum(weights(cols));
end

stats = struct();
stats.area = sum(weights);
stats.normal_flux = normal * weights(:);
stats.patch_areas = patchAreas;
stats.per_face_tags = [];
stats.per_face_areas = [];
stats.per_coarse_triangle_ids = [];
stats.per_coarse_triangle_areas = [];
if isfield(mesh, 'tri_face_tags')
    faceTags = double(mesh.tri_face_tags(:));
    uniqueFaceTags = unique(faceTags, 'stable');
    perFaceAreas = zeros(numel(uniqueFaceTags), 1);
    for k = 1:numel(uniqueFaceTags)
        perFaceAreas(k) = sum(patchAreas(faceTags == uniqueFaceTags(k)));
    end
    stats.per_face_tags = uniqueFaceTags(:);
    stats.per_face_areas = perFaceAreas(:);
end
if isfield(mesh, 'parent_coarse_tri')
    coarseIds = double(mesh.parent_coarse_tri(:));
    uniqueCoarseIds = unique(coarseIds, 'stable');
    perCoarseAreas = zeros(numel(uniqueCoarseIds), 1);
    for k = 1:numel(uniqueCoarseIds)
        perCoarseAreas(k) = sum(patchAreas(coarseIds == uniqueCoarseIds(k)));
    end
    stats.per_coarse_triangle_ids = uniqueCoarseIds(:);
    stats.per_coarse_triangle_areas = perCoarseAreas(:);
end
stats.max_projection_distance = get_meta_scalar(mesh, 'max_projection_distance');
stats.mean_projection_distance = get_meta_scalar(mesh, 'mean_projection_distance');
stats.cad_area = get_meta_scalar(mesh, 'cad_area');
stats.relative_area_error = NaN;
if isfinite(stats.cad_area) && stats.cad_area > 0
    stats.relative_area_error = abs(stats.area - stats.cad_area) / stats.cad_area;
end
end

function value = get_meta_scalar(mesh, field)
value = NaN;
if isfield(mesh, 'meta') && isfield(mesh.meta, field)
    value = double(mesh.meta.(field));
end
end

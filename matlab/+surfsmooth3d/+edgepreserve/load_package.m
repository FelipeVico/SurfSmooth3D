function pkg = load_package(packageDir)
%LOAD_PACKAGE Read and check the fixed CAD, coarse scaffold, and edge source.
packageDir = char(java.io.File(packageDir).getCanonicalPath());
names = {'surface.go3','scaffold.gidmsh','edges_gll.txt','metadata.txt'};
paths = cellfun(@(s) fullfile(packageDir,s),names,'UniformOutput',false);
for k = 1:numel(paths)
    if ~isfile(paths{k}), error('edgepreserve:missingFile', 'Missing %s.',paths{k}); end
end
lines = readlines(paths{4});
meta = struct();
if strtrim(lines(1)) ~= "STEP_MESHER_GEOM_PACKAGE_V1"
    error('edgepreserve:metadata', 'Unsupported geometry package metadata.');
end
for k = 2:numel(lines)
    tokens = regexp(char(lines(k)), '^\s*(\w+)\s+(.+?)\s*$','tokens','once');
    if ~isempty(tokens), meta.(tokens{1}) = tokens{2}; end
end
if ~isfield(meta,'refinement_level') || ~isfield(meta,'surface_order')
    error('edgepreserve:metadata', 'Source order/refinement metadata is required.');
end
sourceLevel = str2double(meta.refinement_level);
sourceOrder = str2double(meta.surface_order);
validateattributes(sourceLevel,{'numeric'},{'scalar','finite','integer','nonnegative'});
validateattributes(sourceOrder,{'numeric'},{'scalar','integer','>=',1,'<=',20});
surface = surfsmooth3d.surfer.load_from_file(paths{1});
[nodes,tris] = surfsmooth3d.edgepreserve.read_gidmsh(paths{2});
if any(~isfinite(nodes(:))) || any(tris(:)<1 | tris(:)>size(nodes,2) | ...
        tris(:) ~= fix(tris(:)))
    error('edgepreserve:scaffold', 'Invalid scaffold coordinates/connectivity.');
end
if any(surface.iptype(:) ~= 1) || any(surface.norders(:) ~= sourceOrder) || ...
        surface.npatches ~= size(tris,2)*4^sourceLevel
    error('edgepreserve:sourceLayout', ...
        'Expected uniform RV patches: source patches = coarse triangles * 4^source refinement.');
end
sourceReport = surfsmooth3d.edgepreserve.validity_report(surface,'fixed CAD source');
if ~sourceReport.is_valid || any(surface.wts <= 0)
    error('edgepreserve:invalidSource', ...
        'CAD source has invalid geometry/normals/weights. Bad nodes: %d.', ...
        sourceReport.stats.bad_node_count);
end
edges = surfsmooth3d.edgepreserve.load_edges(paths{3});
if ~isempty(edges) && isfield(edges,'format_version') && ...
        any([edges.format_version] == 1)
    warning('edgepreserve:legacyEdges', ...
        'V1 edge IDs refer to individual panels, not grouped CAD edges. Prefer a V2 package.');
end
ids = [edges.edge_id];
if numel(unique(ids)) ~= numel(ids)
    error('edgepreserve:edgeIDs', 'CAD edge IDs must be unique.');
end
for k = 1:numel(edges)
    for j = 1:numel(edges(k).panels)
        panel = edges(k).panels(j);
        if any(~isfinite(panel.xyz(:))) || any(~isfinite(panel.w(:))) || any(panel.w < 0)
            error('edgepreserve:edgeQuadrature', 'Invalid GLL edge quadrature.');
        end
    end
end
[~,name] = fileparts(packageDir);
geometryName = name;
if isfield(meta,'source_step'), [~,geometryName] = fileparts(meta.source_step); end
hashes = cellfun(@surfsmooth3d.edgepreserve.file_sha256,paths,'UniformOutput',false);
pkg = struct('package_dir',packageDir,'surface_file',paths{1}, ...
    'scaffold_file',paths{2},'edges_file',paths{3},'metadata_file',paths{4}, ...
    'name',name,'geometry_name',geometryName,'metadata',meta, ...
    'source_paths',{paths},'source_hashes',{hashes},'surface',surface, ...
    'source_order',sourceOrder,'source_level',sourceLevel,'edges',edges, ...
    'coarse_triangle_count',size(tris,2),'scaffold',struct('nodes',nodes,'tris',tris));
end

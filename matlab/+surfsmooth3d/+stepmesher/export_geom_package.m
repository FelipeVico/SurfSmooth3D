function payload = export_geom_package(mesh, packageDir, packageName)
%EXPORT_GEOM_PACKAGE Write a self-contained geometry package folder.
%
% The package contains the fmm3dbie-compatible surface.go3, the flat
% scaffold.gidmsh mesh, one edges_gll.txt file with CAD edges represented
% by piecewise GLL panels, and a small metadata.txt file.

if nargin < 3 || isempty(packageName)
    [~, packageName] = fileparts(char(packageDir));
end

if ~exist(packageDir, 'dir')
    mkdir(packageDir);
end

surfaceFile = fullfile(packageDir, 'surface.go3');
scaffoldFile = fullfile(packageDir, 'scaffold.gidmsh');
edgesFile = fullfile(packageDir, 'edges_gll.txt');
metadataFile = fullfile(packageDir, 'metadata.txt');

surfacePayload = surfsmooth3d.stepmesher.export_go3(mesh, surfaceFile);
scaffoldPayload = surfsmooth3d.stepmesher.export_gidmsh(mesh, scaffoldFile);
edgesPayload = surfsmooth3d.stepmesher.export_edges_gll(mesh, edgesFile);

fid = fopen(metadataFile, 'w');
if fid < 0
    error('stepmesher:io', 'Could not open %s for writing.', metadataFile);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'STEP_MESHER_GEOM_PACKAGE_V1\n');
fprintf(fid, 'name %s\n', packageName);
fprintf(fid, 'source_step %s\n', metadata_string(mesh, 'step_file', ''));
fprintf(fid, 'surface_file surface.go3\n');
fprintf(fid, 'scaffold_file scaffold.gidmsh\n');
fprintf(fid, 'edges_file edges_gll.txt\n');
fprintf(fid, 'surface_order %d\n', surfacePayload.norder);
fprintf(fid, 'edge_order %d\n', edgesPayload.edge_order);
fprintf(fid, 'mesh_fraction %.16e\n', metadata_number(mesh, ...
    'mesh_fraction', NaN));
fprintf(fid, 'refinement_level %d\n', round(metadata_number(mesh, ...
    'refinement_level', -1)));
fprintf(fid, 'num_surface_patches %d\n', surfacePayload.npatches);
fprintf(fid, 'num_cad_edges %d\n', edgesPayload.num_cad_edges);
fprintf(fid, 'num_edge_panels %d\n', edgesPayload.num_panels);
fprintf(fid, 'END_STEP_MESHER_GEOM_PACKAGE\n');

payload = struct();
payload.package_dir = packageDir;
payload.package_name = packageName;
payload.surface_file = surfaceFile;
payload.scaffold_file = scaffoldFile;
payload.edges_file = edgesFile;
payload.metadata_file = metadataFile;
payload.surface = surfacePayload;
payload.scaffold = scaffoldPayload;
payload.edges = edgesPayload;
end

function value = metadata_string(mesh, name, fallback)
value = fallback;
if isfield(mesh, 'meta') && isfield(mesh.meta, name) && ...
        ~isempty(mesh.meta.(name))
    value = char(string(mesh.meta.(name)));
end
end

function value = metadata_number(mesh, name, fallback)
value = fallback;
if isfield(mesh, 'meta') && isfield(mesh.meta, name) && ...
        ~isempty(mesh.meta.(name))
    value = double(mesh.meta.(name));
end
end

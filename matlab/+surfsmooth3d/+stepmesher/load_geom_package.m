function pkg = load_geom_package(packageDir)
%LOAD_GEOM_PACKAGE Load metadata and edges from a geometry package folder.

metadataFile = fullfile(packageDir, 'metadata.txt');
if ~isfile(metadataFile)
    error('stepmesher:fileNotFound', ...
        'Geometry package metadata not found: %s', metadataFile);
end

meta = read_metadata(metadataFile);
pkg = struct();
pkg.package_dir = packageDir;
pkg.metadata_file = metadataFile;
pkg.metadata = meta;
pkg.surface_file = fullfile(packageDir, meta.surface_file);
pkg.scaffold_file = fullfile(packageDir, meta.scaffold_file);
pkg.edges_file = fullfile(packageDir, meta.edges_file);
pkg.edges = surfsmooth3d.stepmesher.load_edges_gll(pkg.edges_file);
if isempty(pkg.edges)
    pkg.edges_format_version = NaN;
else
    pkg.edges_format_version = pkg.edges(1).format_version;
end

if ~isempty(which('surfsmooth3d.surfer.load_from_file')) && isfile(pkg.surface_file)
    pkg.surface = surfsmooth3d.surfer.load_from_file(pkg.surface_file);
end
end

function meta = read_metadata(metadataFile)
fid = fopen(metadataFile, 'r');
if fid < 0
    error('stepmesher:io', 'Could not open %s for reading.', metadataFile);
end
cleanup = onCleanup(@() fclose(fid));

line = next_nonempty_line(fid);
if ~strcmp(line, 'STEP_MESHER_GEOM_PACKAGE_V1')
    error('stepmesher:badPackage', ...
        'Unsupported geometry package metadata marker.');
end

meta = struct();
while true
    line = next_nonempty_line(fid);
    if strcmp(line, 'END_STEP_MESHER_GEOM_PACKAGE')
        break
    end
    parts = strsplit(line);
    key = matlab.lang.makeValidName(parts{1});
    value = strjoin(parts(2:end), ' ');
    numericValue = str2double(value);
    if ~isnan(numericValue)
        meta.(key) = numericValue;
    else
        meta.(key) = char(value);
    end
end
end

function line = next_nonempty_line(fid)
while true
    line = fgetl(fid);
    if ~ischar(line)
        error('stepmesher:badPackage', 'Unexpected end of file.');
    end
    line = strtrim(line);
    if ~isempty(line)
        return
    end
end
end

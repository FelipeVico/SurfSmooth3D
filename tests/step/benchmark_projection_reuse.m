function capture = benchmark_projection_reuse(packageRoot, outputFile, selectedCases, settings)
%BENCHMARK_PROJECTION_REUSE Capture full meshes from one isolated build.
% Run baseline and candidate in SEPARATE MATLAB processes. packageRoot is the
% CMake build/install root containing matlab/, not the source repository.
% selectedCases is a cell/string list of STEP stems, 'matrix' (default), or
% 'all'. settings accepts orders (default 4), hairy_orders (default 4),
% repetitions (default 1), step_directory, and label.
% The native mesh call is timed. Diagnostic builds include request/BRep writes
% in that call; their timing CSV separates this diagnostic_capture cost.
% srcvals, area/flux validation, MAT checkpoints, and JSON summaries are
% outside that timer; no plots or geometry exports run.

if nargin < 3 || isempty(selectedCases), selectedCases = 'matrix'; end
if nargin < 4 || isempty(settings), settings = struct(); end
sourceRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
settings = defaults(settings, sourceRoot);
packageRoot = canonical(packageRoot);
matlabRoot = fullfile(packageRoot, 'matlab');
assert(isfolder(matlabRoot), 'stepmesher:benchmarkPackage', ...
    'Build root must contain matlab/: %s', packageRoot);
outputFile = canonical(outputFile);
outputDir = fileparts(outputFile);
if ~isfolder(outputDir), mkdir(outputDir); end

originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
restoredefaultpath;
addpath(matlabRoot, '-begin');
resolved = verify_package(matlabRoot);
names = select_cases(selectedCases, settings.step_directory);

capture = struct('schema_version', 1, 'label', settings.label, ...
    'package_root', packageRoot, 'resolved_functions', resolved, ...
    'matlab_version', version, 'platform', computer, ...
    'created_utc', char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
    'settings', settings, 'cases', {{}});
fprintf('PROJECTION_BENCHMARK_PACKAGE %s\n', packageRoot);
fprintf('PROJECTION_BENCHMARK_MEX %s\n', resolved.expected_mex);
save(outputFile, 'capture', '-v7');

for i = 1:numel(names)
    name = names{i};
    if is_hairy(name), orders = settings.hairy_orders;
    else, orders = settings.orders; end
    for order = orders(:).'
        [options, choice] = case_options(name, order);
        entry = struct('id', sprintf('%s_p%d_%s', name, order, choice), ...
            'name', name, 'options', options, 'profile_choice', choice, ...
            'runs', {{}}, 'passed', true, 'error_identifier', '', ...
            'error_message', '');
        stepFile = fullfile(settings.step_directory, [name, '.step']);
        assert(isfile(stepFile), 'stepmesher:benchmarkFixture', ...
            'Missing STEP fixture: %s', stepFile);
        for repeat = 1:settings.repetitions
            fprintf('PROJECTION_BENCHMARK_BEGIN %s repetition=%d\n', entry.id, repeat);
            try
                capturePrefix = fullfile(outputDir, sprintf('%s_p%d', name, order));
                if settings.repetitions > 1
                    capturePrefix = sprintf('%s_r%d', capturePrefix, repeat);
                end
                [mesh, meshSeconds] = timed_mesh(stepFile, options, capturePrefix);
                loadedMex = verify_loaded_mex(resolved.expected_mex);
                timer = tic;
                [srcvals, norders, iptype] = surfsmooth3d.stepmesher.to_srcvals(mesh);
                srcvalsSeconds = toc(timer);
                timer = tic;
                stats = surfsmooth3d.stepmesher.validate_area_flux(mesh);
                validationSeconds = toc(timer);
                run = struct('mesh', mesh, 'srcvals', srcvals, ...
                    'norders', norders, 'iptype', iptype, 'stats', stats, ...
                    'loaded_mex', loadedMex, 'mesh_seconds', meshSeconds, ...
                    'srcvals_seconds', srcvalsSeconds, ...
                    'validation_seconds', validationSeconds, ...
                    'native_capture', native_capture(capturePrefix));
                entry.runs{end+1} = run;
                fprintf(['PROJECTION_BENCHMARK_RESULT %s repetition=%d ', ...
                    'seconds=%.9g patches=%d nodes=%d failures=%d ', ...
                    'area_relative_error=%.9g\n'], entry.id, repeat, ...
                    meshSeconds, mesh.npatches, size(mesh.xyz, 2), ...
                    mesh.meta.projection_failures, stats.relative_area_error);
            catch error
                entry.passed = false;
                entry.error_identifier = error.identifier;
                entry.error_message = error.message;
                fprintf(2, 'PROJECTION_BENCHMARK_FAILED %s: %s\n', entry.id, error.message);
                break;
            end
            capture.cases{end+1} = entry;
            checkpoint(outputFile, capture);
            capture.cases(end) = [];
        end
        capture.cases{end+1} = entry;
        checkpoint(outputFile, capture);
    end
end
assert(all(cellfun(@(entry) entry.passed, capture.cases)), ...
    'stepmesher:benchmarkFailed', 'A benchmark case failed; inspect %s.', outputFile);
fprintf('PROJECTION_BENCHMARK_COMPLETE %s\n', outputFile);
end

function settings = defaults(settings, root)
values = struct('orders', 4, 'hairy_orders', 4, 'repetitions', 1, ...
    'step_directory', fullfile(root, 'step_mesher', 'examples', 'step_files'), ...
    'label', '');
unknown = setdiff(fieldnames(settings), fieldnames(values));
assert(isempty(unknown), 'stepmesher:benchmarkSettings', ...
    'Unknown settings: %s', strjoin(unknown, ', '));
fields = fieldnames(values);
for i = 1:numel(fields)
    if ~isfield(settings, fields{i}), settings.(fields{i}) = values.(fields{i}); end
end
for name = {'orders', 'hairy_orders'}
    orders = settings.(name{1});
    assert(isnumeric(orders) && isvector(orders) && ~isempty(orders) && ...
        all(isfinite(orders)) && all(orders >= 1) && all(orders == fix(orders)), ...
        'stepmesher:benchmarkOrder', 'Orders must be positive integers.');
    settings.(name{1}) = unique(orders(:).', 'stable');
end
assert(isscalar(settings.repetitions) && isfinite(settings.repetitions) && ...
    settings.repetitions >= 1 && settings.repetitions == fix(settings.repetitions), ...
    'stepmesher:benchmarkRepetitions', 'Repetitions must be a positive integer.');
settings.step_directory = canonical(settings.step_directory);
settings.label = char(settings.label);
end

function names = select_cases(selection, stepDirectory)
if ischar(selection) || (isstring(selection) && isscalar(selection))
    selection = char(selection);
    if strcmp(selection, 'matrix')
        names = {'test_geometry_1', 'test_curved_uv', 'spheres_intersect', ...
            'Multi_res', 'UV_curved_bezier', 'piecewise_bezier_cap', ...
            'torus_ellipse', 'wobbly_hairy_torus_10_v2_flat_AP214'};
        return;
    elseif strcmp(selection, 'all')
        files = dir(fullfile(stepDirectory, '*.step'));
        names = cellfun(@(name) erase(name, '.step'), {files.name}, 'UniformOutput', false);
        names = sort(names);
        return;
    else
        selection = {selection};
    end
end
if isstring(selection), selection = cellstr(selection); end
assert(iscell(selection) && ~isempty(selection), ...
    'stepmesher:benchmarkCases', 'Select at least one STEP fixture.');
names = cellfun(@char, selection, 'UniformOutput', false);
names = cellfun(@(name) erase(name, '.step'), names, 'UniformOutput', false);
assert(all(cellfun(@(name) isempty(regexp(name, '[/\\]', 'once')), names)), ...
    'stepmesher:benchmarkCases', 'Case names must be STEP filename stems.');
assert(numel(unique(names)) == numel(names), ...
    'stepmesher:benchmarkCases', 'Duplicate STEP cases are not allowed.');
end

function [options, choice] = case_options(name, order)
options = struct('order', order, 'mesh_fraction', 0.10, ...
    'refinement_level', 0, 'edge_gll_order', order);
if is_hairy(name), choice = 'curvature_balanced_local';
else, choice = 'rigid'; end
[options, ~] = surfsmooth3d.stepmesher.apply_mesh_profile_choice(options, choice);
% These match the running example, including its explicit import settings.
options.occt_precision = 1e-6;
options.occt_maxprecision = 1e-6;
options.sameparameter = true;
options.surfacecurve_mode = '3d_preferred';
options.verbose = false;
end

function yes = is_hairy(name)
yes = strcmp(name, 'wobbly_hairy_torus_10_v2_flat_AP214');
end

function [mesh, seconds] = timed_mesh(stepFile, options, capturePrefix)
previous = getenv('SS3D_PROJECTION_CAPTURE_PREFIX');
cleanup = onCleanup(@() setenv('SS3D_PROJECTION_CAPTURE_PREFIX', previous));
setenv('SS3D_PROJECTION_CAPTURE_PREFIX', capturePrefix);
timer = tic;
mesh = surfsmooth3d.stepmesher.mesh_step(stepFile, options);
seconds = toc(timer);
end

function value = native_capture(prefix)
value = struct('prefix', prefix, 'request_file', [prefix, '.request.txt'], ...
    'map_file', [prefix, '.map.txt'], 'brep_file', [prefix, '.brep'], ...
    'timings_file', [prefix, '.timings.csv'], 'timings_csv', '');
if isfile(value.timings_file), value.timings_csv = fileread(value.timings_file); end
end

function resolved = verify_package(matlabRoot)
names = {'surfsmooth3d.stepmesher.mesh_step', ...
    'surfsmooth3d.stepmesher.to_srcvals', ...
    'surfsmooth3d.stepmesher.validate_area_flux', ...
    'surfsmooth3d.stepmesher.apply_mesh_profile_choice', ...
    'surfsmooth3d.internal.koorn.rv_nodes', ...
    'surfsmooth3d.internal.koorn.rv_weights', ...
    'surfsmooth3d.internal.koorn.ders', ...
    'surfsmooth3d.internal.koorn.vals2coefs'};
resolved = struct('names', {names}, 'paths', {cell(size(names))}, ...
    'expected_mex', fullfile(matlabRoot, '+surfsmooth3d', '+stepmesher', ...
        'private', ['step_mesher_mex.', mexext]));
for i = 1:numel(names)
    found = which(names{i});
    assert(~isempty(found), 'stepmesher:benchmarkResolution', 'Missing %s.', names{i});
    found = canonical(found);
    assert(startsWith(found, [matlabRoot, filesep]), ...
        'stepmesher:benchmarkShadow', '%s resolved outside build package: %s', names{i}, found);
    resolved.paths{i} = found;
end
assert(isfile(resolved.expected_mex), 'stepmesher:benchmarkMex', ...
    'Missing expected private MEX: %s', resolved.expected_mex);
resolved.expected_mex = canonical(resolved.expected_mex);
[~, mexFiles] = inmem('-completenames');
for i = 1:numel(mexFiles)
    [~, stem] = fileparts(mexFiles{i});
    if strcmp(stem, 'step_mesher_mex')
        assert(strcmp(canonical(mexFiles{i}), resolved.expected_mex), ...
            'stepmesher:benchmarkMexLoaded', ...
            'Another STEP MEX is already loaded; use a separate MATLAB process.');
    end
end
end

function loaded = verify_loaded_mex(expected)
[~, mexFiles] = inmem('-completenames');
loaded = '';
for i = 1:numel(mexFiles)
    [~, stem] = fileparts(mexFiles{i});
    if strcmp(stem, 'step_mesher_mex')
        found = canonical(mexFiles{i});
        assert(strcmp(found, expected), 'stepmesher:benchmarkMexShadow', ...
            'STEP MEX loaded from the wrong build: %s', found);
        loaded = found;
    end
end
assert(~isempty(loaded), 'stepmesher:benchmarkMexLoaded', ...
    'Expected STEP MEX was not found among loaded MEX files.');
end

function checkpoint(outputFile, capture)
save(outputFile, 'capture', '-v7');
summary = rmfield(capture, 'cases');
summary.cases = cell(size(capture.cases));
for i = 1:numel(capture.cases)
    entry = capture.cases{i};
    rows = cell(size(entry.runs));
    for j = 1:numel(entry.runs)
        run = entry.runs{j};
        rows{j} = struct('mesh_seconds', run.mesh_seconds, ...
            'srcvals_seconds', run.srcvals_seconds, ...
            'validation_seconds', run.validation_seconds, ...
            'native_capture', run.native_capture, ...
            'loaded_mex', run.loaded_mex, 'patches', run.mesh.npatches, ...
            'nodes', size(run.mesh.xyz, 2), ...
            'projection_failures', run.mesh.meta.projection_failures, ...
            'area', run.stats.area, 'relative_area_error', run.stats.relative_area_error, ...
            'normal_flux', run.stats.normal_flux, ...
            'max_projection_distance', run.stats.max_projection_distance, ...
            'mean_projection_distance', run.stats.mean_projection_distance);
    end
    entry.runs = rows;
    summary.cases{i} = entry;
end
[folder, stem] = fileparts(outputFile);
fid = fopen(fullfile(folder, [stem, '.json']), 'w');
assert(fid >= 0, 'stepmesher:benchmarkOutput', 'Could not write benchmark summary.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(summary, PrettyPrint=true));
end

function value = canonical(value)
value = char(java.io.File(char(value)).getCanonicalPath());
end

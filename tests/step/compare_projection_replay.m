function report = compare_projection_replay(candidateFile, replayReportFile, outputFile)
%COMPARE_PROJECTION_REPLAY Check a continuous MATLAB mesh against full replay.
% Reads existing captures and replay files only; does not invoke the mesher.
% The complete replay must have passed strict baseline/candidate byte checks.
% Its XYZ/distance rows are reconstructed in original order and compared
% bitwise to the actual MATLAB mesh, whose face cache runs continuously.
% Projected UV is absent from the MATLAB mesh; replay byte checks cover UV.
if nargin < 3, outputFile = ''; end
replay = jsondecode(fileread(replayReportFile));
assert(replay.completed && replay.passed && ...
    replay.verified_nodes == replay.source_nodes && isempty(replay.errors), ...
    'stepmesher:replayIncomplete', 'Full replay is incomplete or failed.');
assert(strcmp(file_sha256(replay.request), replay.request_sha256), ...
    'stepmesher:replayProvenance', 'Original replay request hash changed.');
request = read_request(replay.request);
assert(request.count == replay.source_nodes, 'stepmesher:replayCount', ...
    'Original request count differs from the replay report.');
count = request.count;
expectedXYZ = zeros(3, count);
expectedDistance = zeros(1, count);
fallbacks = 0;
rows = replay.batches;
% jsondecode uses cells when resumed batches contain additional fields.
if isstruct(rows), rows = num2cell(rows); end
assert(iscell(rows) && all(cellfun(@isstruct, rows)), ...
    'stepmesher:replayFormat', 'Replay batches must be structures.');
indexes = cellfun(@(batch) batch.index, rows);
[~, order] = sort(indexes);
rows = rows(order);
assert(isequal(reshape(indexes(order), 1, []), 0:numel(rows)-1), ...
    'stepmesher:replayCoverage', 'Replay batch indexes are incomplete.');
next = 0;
replayDir = fileparts(replayReportFile);
for i = 1:numel(rows)
    batch = rows{i};
    assert(batch.exact && batch.begin == next && batch.end > batch.begin && ...
        batch.end <= count && batch.nodes == batch.end - batch.begin, ...
        'stepmesher:replayCoverage', 'Replay ranges are incomplete or overlapping.');
    batchDir = fullfile(replayDir, sprintf('batch-%04d', batch.index));
    baselineFile = fullfile(batchDir, 'baseline.txt');
    candidateReplayFile = fullfile(batchDir, 'candidate.txt');
    assert(strcmp(file_sha256(baselineFile), batch.baseline_sha256) && ...
        strcmp(file_sha256(candidateReplayFile), batch.candidate_sha256) && ...
        isequal(file_bytes(baselineFile), file_bytes(candidateReplayFile)), ...
        'stepmesher:replayProvenance', 'Replay response bytes or hashes differ.');
    indices = batch.begin+1:batch.end;
    chunkRequest = read_request(fullfile(batchDir, 'request.txt'));
    assert(strcmp(chunkRequest.step_file, request.step_file) && ...
        same_bits(chunkRequest.faces, request.faces) && ...
        chunkRequest.count == batch.nodes && ...
        same_bits(chunkRequest.nodes, request.nodes(:, indices)), ...
        'stepmesher:replayOrder', 'Batch request differs from original point or face order.');
    [values, chunkFallbacks] = read_response(baselineFile, batch.nodes);
    assert(chunkFallbacks == batch.fallbacks, 'stepmesher:replayFallback', ...
        'Batch fallback count differs from the replay report.');
    expectedXYZ(:, indices) = values(1:3, :);
    expectedDistance(indices) = values(4, :);
    fallbacks = fallbacks + chunkFallbacks;
    next = batch.end;
end
assert(next == count, 'stepmesher:replayCoverage', ...
    'Replay batches do not cover the whole original request.');

loaded = load(candidateFile, 'capture');
assert(isfield(loaded, 'capture') && ~isempty(loaded.capture.cases), ...
    'stepmesher:replayCapture', 'Candidate capture has no complete mesh.');
capture = loaded.capture;
report = struct('candidate_file', candidateFile, ...
    'candidate_sha256', file_sha256(candidateFile), ...
    'replay_report_file', replayReportFile, ...
    'replay_report_sha256', file_sha256(replayReportFile), ...
    'request_file', replay.request, ...
    'request_sha256', replay.request_sha256, 'nodes', count, ...
    'batches', numel(rows), 'fallbacks', fallbacks, ...
    'replay_xyz_distance_uv_byte_exact', true, ...
    'projected_uv_in_matlab_mesh', false, 'passed', true, 'runs', {{}});
for i = 1:numel(capture.cases)
    entry = capture.cases{i};
    for j = 1:numel(entry.runs)
        run = entry.runs{j};
        % The instrumentation writes this run's actual preprojection inputs.
        % Hash equality links a complete mesh to the same replay input batch.
        if ~isfield(run, 'native_capture') || ...
                ~isfile(run.native_capture.request_file) || ...
                ~strcmp(file_sha256(run.native_capture.request_file), replay.request_sha256)
            continue;
        end
        assert(entry.passed && isfield(run, 'mesh'), 'stepmesher:replayCapture', ...
            'Matching candidate run did not complete successfully.');
        mesh = run.mesh;
        assert(isfield(mesh, 'proj_dist') && isfield(mesh, 'tri_face_tags') && ...
            isfield(mesh, 'nodes_per_patch') && isfield(mesh, 'meta') && ...
            isfield(mesh.meta, 'projection_failures'), ...
            'stepmesher:replayCapture', 'Candidate mesh lacks projection or face fields.');
        tags = repelem(reshape(mesh.tri_face_tags, 1, []), mesh.nodes_per_patch);
        xyzExact = same_bits(mesh.xyz, expectedXYZ);
        distanceExact = same_bits(reshape(mesh.proj_dist, 1, []), expectedDistance);
        faceOrderExact = same_bits(tags, request.nodes(1, :));
        fallbackExact = mesh.meta.projection_failures == fallbacks;
        shapeExact = size(mesh.xyz, 1) == 3 && size(mesh.xyz, 2) == count && ...
            numel(mesh.proj_dist) == count && numel(tags) == count && ...
            mesh.npatches == numel(mesh.tri_face_tags);
        checks = struct('id', entry.id, 'repetition', j, ...
            'xyz_bitwise_exact', xyzExact, 'distance_bitwise_exact', distanceExact, ...
            'face_order_bitwise_exact', faceOrderExact, ...
            'fallback_count_exact', fallbackExact, 'dimensions_exact', shapeExact, ...
            'passed', xyzExact && distanceExact && faceOrderExact && fallbackExact && shapeExact);
        report.runs{end+1} = checks;
        report.passed = report.passed && checks.passed;
        fprintf('PROJECTION_REPLAY_CONTINUOUS %s repetition=%d nodes=%d exact=%d\n', ...
            entry.id, j, count, checks.passed);
    end
end
assert(~isempty(report.runs), 'stepmesher:replayCapture', ...
    'No complete MATLAB mesh is linked to the original replay request.');
if ~isempty(outputFile)
    fid = fopen(outputFile, 'w');
    assert(fid >= 0, 'stepmesher:replayOutput', 'Could not write comparison report.');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', jsonencode(report, PrettyPrint=true));
end
assert(report.passed, 'stepmesher:replayChanged', ...
    'Continuous MATLAB mesh differs from baseline replay; inspect the comparison report.');
end

function request = read_request(filename)
fid = fopen(filename, 'r');
assert(fid >= 0, 'stepmesher:replayRead', 'Could not read %s.', filename);
cleanup = onCleanup(@() fclose(fid));
assert(strcmp(fgetl(fid), 'STEP_MESHER_OCC_PROJECT_V1'), ...
    'stepmesher:replayFormat', 'Unexpected projection request format.');
request.step_file = fgetl(fid);
faceCount = read_count(fid);
request.faces = zeros(8, faceCount);
for i = 1:faceCount
    values = sscanf(fgetl(fid), '%f');
    assert(numel(values) == 8 && all(isfinite(values)), ...
        'stepmesher:replayFormat', 'Invalid face descriptor.');
    request.faces(:, i) = values;
end
request.count = read_count(fid);
request.nodes = read_numbers(fid, 4, request.count);
assert(all(request.nodes(1, :) > 0 & request.nodes(1, :) == fix(request.nodes(1, :))), ...
    'stepmesher:replayFormat', 'Invalid original face tags.');
end

function [values, fallbacks] = read_response(filename, expectedCount)
fid = fopen(filename, 'r');
assert(fid >= 0, 'stepmesher:replayRead', 'Could not read %s.', filename);
cleanup = onCleanup(@() fclose(fid));
header = sscanf(fgetl(fid), 'OK %d %d');
assert(numel(header) == 2 && header(1) == expectedCount && ...
    header(2) >= 0 && header(2) <= expectedCount, ...
    'stepmesher:replayFormat', 'Unexpected projection response header.');
fallbacks = header(2);
values = read_numbers(fid, 6, expectedCount);
end

function count = read_count(fid)
value = sscanf(fgetl(fid), '%f');
assert(isscalar(value) && isfinite(value) && value >= 0 && value == fix(value), ...
    'stepmesher:replayFormat', 'Invalid request count.');
count = value;
end

function values = read_numbers(fid, width, count)
[values, readCount] = fscanf(fid, '%f', [width, count]);
assert(readCount == width*count && all(isfinite(values(:))) && ...
    isempty(strtrim(fread(fid, Inf, '*char').')), ...
    'stepmesher:replayFormat', 'Incomplete, nonfinite, or extra projection values.');
end

function yes = same_bits(a, b)
yes = isa(a, 'double') && isa(b, 'double') && isequal(size(a), size(b)) && ...
    isequal(typecast(a(:), 'uint8'), typecast(b(:), 'uint8'));
end

function bytes = file_bytes(filename)
fid = fopen(filename, 'rb');
assert(fid >= 0, 'stepmesher:replayRead', 'Could not read %s.', filename);
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
end

function value = file_sha256(filename)
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(typecast(file_bytes(filename), 'int8'));
bytes = typecast(digest.digest(), 'uint8');
value = lower(reshape(dec2hex(bytes, 2).', 1, []));
end

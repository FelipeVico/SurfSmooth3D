function report = compare_projection_reuse(baselineFile, candidateFile, outputFile, requireExact)
%COMPARE_PROJECTION_REUSE Compare full benchmark captures without tolerances.
% Numeric isequaln equality and byte equality are recorded separately: NaN
% payloads or signed zero can differ in bytes while MATLAB values compare
% equal. Only step_file and occ_info_exe paths are excluded from mesh payloads.
% The default acceptance gate requires both numeric and byte equality.
% Timing is reported independently and never weakens numerical acceptance.
if nargin < 3, outputFile = ''; end
if nargin < 4, requireExact = true; end
baseline = load(baselineFile, 'capture'); baseline = baseline.capture;
candidate = load(candidateFile, 'capture'); candidate = candidate.capture;
assert(baseline.schema_version == candidate.schema_version, ...
    'stepmesher:benchmarkSchema', 'Capture schema versions differ.');
baseIds = cellfun(@(entry) entry.id, baseline.cases, 'UniformOutput', false);
candidateIds = cellfun(@(entry) entry.id, candidate.cases, 'UniformOutput', false);
assert(numel(unique(baseIds)) == numel(baseIds) && ...
    numel(unique(candidateIds)) == numel(candidateIds) && ...
    isequal(sort(baseIds), sort(candidateIds)), ...
    'stepmesher:benchmarkCases', 'Baseline and candidate case sets differ.');
report = struct('baseline_file', baselineFile, 'candidate_file', candidateFile, ...
    'passed', true, 'bitwise_passed', true, 'cases', {{}});
for i = 1:numel(baseIds)
    expected = baseline.cases{i};
    actual = candidate.cases{find(strcmp(candidateIds, baseIds{i}), 1)};
    comparisons = {};
    comparisons = compare_value(expected.options, actual.options, 'options', comparisons);
    comparisons = compare_value(expected.passed, actual.passed, 'capture.passed', comparisons);
    if ~expected.passed || ~actual.passed || isempty(expected.runs) || isempty(actual.runs)
        error('stepmesher:benchmarkIncomplete', 'Case %s has no successful complete capture.', expected.id);
    end
    % Every candidate repetition and every baseline repetition must reproduce
    % the first baseline payload. Timings and build paths are not payloads.
    reference = payload(expected.runs{1});
    for j = 1:numel(expected.runs)
        comparisons = compare_value(reference, payload(expected.runs{j}), ...
            sprintf('baseline_repeat_%d', j), comparisons);
    end
    for j = 1:numel(actual.runs)
        comparisons = compare_value(reference, payload(actual.runs{j}), ...
            sprintf('candidate_repeat_%d', j), comparisons);
    end
    exact = all(cellfun(@(entry) entry.isequaln, comparisons));
    bitwise = all(cellfun(@(entry) entry.bitwise, comparisons));
    baseSeconds = cellfun(@(run) run.mesh_seconds, expected.runs);
    candidateSeconds = cellfun(@(run) run.mesh_seconds, actual.runs);
    entry = struct('id', expected.id, 'passed', exact, 'bitwise_passed', bitwise, ...
        'baseline_seconds', baseSeconds, 'candidate_seconds', candidateSeconds, ...
        'baseline_median_seconds', median(baseSeconds), ...
        'candidate_median_seconds', median(candidateSeconds), ...
        'speedup', median(baseSeconds) / median(candidateSeconds), ...
        'comparisons', {comparisons});
    report.cases{end+1} = entry;
    report.passed = report.passed && exact;
    report.bitwise_passed = report.bitwise_passed && bitwise;
    fprintf('PROJECTION_REUSE_COMPARISON %s exact=%d bitwise=%d speedup=%.6g\n', ...
        expected.id, exact, bitwise, entry.speedup);
    changed = comparisons(~cellfun(@(value) value.isequaln && value.bitwise, comparisons));
    for j = 1:min(10, numel(changed))
        fprintf(2, '  DIFFERENCE %s exact=%d bitwise=%d absolute=%.17g scaled=%.17g\n', ...
            changed{j}.label, changed{j}.isequaln, changed{j}.bitwise, ...
            changed{j}.maximum_absolute_difference, ...
            changed{j}.maximum_scaled_difference);
    end
end
if ~isempty(outputFile)
    folder = fileparts(outputFile);
    if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
    fid = fopen(outputFile, 'w');
    assert(fid >= 0, 'stepmesher:benchmarkOutput', 'Could not write comparison report.');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', jsonencode(report, PrettyPrint=true));
end
if requireExact
    assert(report.passed && report.bitwise_passed, 'stepmesher:benchmarkChanged', ...
        'Candidate or repeated baseline output values or bytes differ; inspect the comparison report.');
end
end

function value = payload(run)
value = struct('mesh', run.mesh, 'srcvals', run.srcvals, ...
    'norders', run.norders, 'iptype', run.iptype, 'stats', run.stats);
end

function comparisons = compare_value(expected, actual, label, comparisons)
if isstruct(expected) && isstruct(actual) && isequal(size(expected), size(actual)) && ...
        isequal(sort(fieldnames(expected)), sort(fieldnames(actual)))
    fields = fieldnames(expected);
    for k = 1:numel(expected)
        for i = 1:numel(fields)
            name = fields{i};
            if any(strcmp(name, {'step_file', 'occ_info_exe'})), continue; end
            next = [label, '.', name];
            if numel(expected) > 1, next = sprintf('%s(%d)', next, k); end
            comparisons = compare_value(expected(k).(name), actual(k).(name), next, comparisons);
        end
    end
    return;
end
if iscell(expected) && iscell(actual) && isequal(size(expected), size(actual))
    for i = 1:numel(expected)
        comparisons = compare_value(expected{i}, actual{i}, ...
            sprintf('%s{%d}', label, i), comparisons);
    end
    return;
end
entry = struct('label', label, 'isequaln', isequaln(expected, actual), ...
    'bitwise', false, 'same_class', strcmp(class(expected), class(actual)), ...
    'same_size', isequal(size(expected), size(actual)), ...
    'same_finite_pattern', true, 'maximum_absolute_difference', 0, ...
    'maximum_scaled_difference', 0);
entry.bitwise = entry.isequaln && entry.same_class && entry.same_size;
if (isnumeric(expected) || islogical(expected)) && ...
        (isnumeric(actual) || islogical(actual)) && entry.same_class && entry.same_size
    entry.bitwise = same_bytes(expected, actual);
    entry.same_finite_pattern = isequal(isfinite(expected), isfinite(actual)) && ...
        isequal(isnan(expected), isnan(actual));
    use = isfinite(expected) & isfinite(actual);
    if any(use(:))
        delta = max(abs(double(expected(use)) - double(actual(use))));
        scale = max(1, max(abs(double(expected(use)))));
        entry.maximum_absolute_difference = delta;
        entry.maximum_scaled_difference = delta / scale;
    end
    if ~entry.same_finite_pattern || ~isequaln(expected(~use), actual(~use))
        entry.maximum_absolute_difference = Inf;
        entry.maximum_scaled_difference = Inf;
    end
elseif ~entry.isequaln
    entry.maximum_absolute_difference = Inf;
    entry.maximum_scaled_difference = Inf;
end
comparisons{end+1} = entry;
end

function yes = same_bytes(expected, actual)
if islogical(expected)
    yes = isequal(expected, actual);
elseif issparse(expected) || issparse(actual)
    yes = isequaln(expected, actual);
else
    yes = isequal(typecast(real(expected(:)), 'uint8'), typecast(real(actual(:)), 'uint8'));
    if ~isreal(expected) || ~isreal(actual)
        yes = yes && isequal(typecast(imag(expected(:)), 'uint8'), ...
            typecast(imag(actual(:)), 'uint8'));
    end
end
end

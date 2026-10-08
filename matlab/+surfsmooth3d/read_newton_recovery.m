function info = read_newton_recovery(filename,enabled,preserve)
%READ_NEWTON_RECOVERY Read tagged NEWTON_RECOVERY_V1 stage/refinement records.
% Missing reports mean no targets were deferred. Optional preserve=true copies
% a present report out of a temporary solver directory before it is removed.
if nargin < 2, enabled = true; end
if nargin < 3, preserve = false; end
info = struct('enabled',logical(enabled),'deferred',0,'recovered',0, ...
    'unresolved',0,'events',struct([]),'report_file','','report_available',false, ...
    'validation_failed',false);
if ~isfile(filename), return; end
fid = fopen(filename,'r');
if fid < 0
    error('MULTISCALE_MESHER:RecoveryDiagnosticFile','Cannot read recovery report %s.',filename);
end
cleanup = onCleanup(@() fclose(fid));
header = fgetl(fid);
if ~ischar(header) || ~strcmp(strtrim(header),'NEWTON_RECOVERY_V1')
    error('MULTISCALE_MESHER:RecoveryDiagnosticFile','Unsupported recovery report %s.',filename);
end
event = []; target = [];
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line), continue; end
    [tag,value] = strtok(line); value = strtrim(value);
    switch tag
        case 'BEGIN'
            numbers = values(value,3);
            event = struct('stage',numbers(1),'refinement',numbers(2), ...
                'deferred',numbers(3),'recovered',0,'unresolved',0,'targets',struct([]));
        case 'TARGET'
            numbers = values(value,3);
            target = struct('index',numbers(1),'reason',numbers(2),'iteration',numbers(3), ...
                'outcome','','base_point',nan(1,3),'normal',nan(1,3), ...
                'initial_point',nan(1,3),'last_point',nan(1,3), ...
                'initial_height',NaN,'last_height',NaN, ...
                'rejected_height',NaN,'rejected_point',nan(1,3),'last_residual',NaN, ...
                'last_gradient',nan(1,3),'sigma',NaN,'interval',nan(1,2), ...
                'root',NaN,'residual',NaN,'slope',NaN,'position_tolerance',NaN, ...
                'counts',zeros(1,3),'brackets',zeros(0,4), ...
                'candidates',zeros(0,4),'candidate_outcomes',{{}});
        case 'OUTCOME'
            target.outcome = value;
        case {'INITIAL_HEIGHT','LAST_HEIGHT','REJECTED_HEIGHT','LAST_RESIDUAL', ...
                'SIGMA','ROOT','RESIDUAL','SLOPE','POSITION_TOLERANCE'}
            target.(lower(tag)) = values(value,1);
        case {'BASE_POINT','NORMAL','INITIAL_POINT','LAST_POINT', ...
                'REJECTED_POINT','LAST_GRADIENT','COUNTS'}
            target.(lower(tag)) = values(value,3);
        case 'INTERVAL'
            target.interval = values(value,2);
        case 'BRACKET'
            target.brackets(end+1,:) = values(value,4);
        case 'CANDIDATE'
            target.candidates(end+1,:) = values(value,4);
        case 'CANDIDATE_OUTCOME'
            target.candidate_outcomes{end+1} = value;
        case 'END_TARGET'
            if isempty(event) || isempty(target), malformed(filename); end
            if isempty(event.targets)
                event.targets = target;
            else
                event.targets(end+1) = target;
            end
            target = [];
        case 'END'
            numbers = values(value,2);
            if isempty(event) || ~isempty(target) || numel(event.targets)~=event.deferred
                malformed(filename);
            end
            event.recovered = numbers(1); event.unresolved = numbers(2);
            if event.recovered+event.unresolved~=event.deferred, malformed(filename); end
            if isempty(info.events)
                info.events = event;
            else
                info.events(end+1) = event;
            end
            info.deferred = info.deferred+event.deferred;
            info.recovered = info.recovered+event.recovered;
            info.unresolved = info.unresolved+event.unresolved;
            event = [];
        case 'VALIDATION_FAILED'
            info.validation_failed = true;
        otherwise
            % Future tagged detail can be added without changing this reader.
    end
end
if ~isempty(event) || ~isempty(target), malformed(filename); end
info.report_file = filename;
info.report_available = true;
if preserve
    destination = [tempname '_newton_recovery.txt'];
    [ok,message] = copyfile(filename,destination);
    if ~ok
        error('MULTISCALE_MESHER:RecoveryDiagnosticFile', ...
            'Could not preserve recovery diagnostic %s: %s',filename,message);
    end
    info.report_file = destination;
end
end

function x = values(text,count)
text = regexprep(text,'[dD]([+-]?\d+)','e$1');
text = strrep(text,'Infinity','Inf');
parts = regexp(strtrim(text),'\s+','split');
if numel(parts)~=count
    error('MULTISCALE_MESHER:RecoveryDiagnosticFile','Malformed recovery numeric record.');
end
x = cellfun(@str2double,parts);
end

function malformed(filename)
error('MULTISCALE_MESHER:RecoveryDiagnosticFile','Incomplete recovery report %s.',filename);
end

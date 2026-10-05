function report = read_newton_failure(filename)
%READ_NEWTON_FAILURE Read a radius-guard or ordinary Newton failure snapshot.
fid = fopen(filename,'r');
if fid < 0, error('edgepreserve:newtonDiagnostic','Cannot open diagnostic: %s',filename); end
cleanup = onCleanup(@() fclose(fid));
version = strtrim(fgetl(fid));
if ~ismember(version,{'NEWTON_RADIUS_FAILURE_V1','NEWTON_PROJECTION_FAILURE_V2'})
    error('edgepreserve:newtonDiagnostic','Unsupported Newton diagnostic: %s',filename);
end
report = struct('filename',filename,'failure_kind','radius', ...
    'failure_reason','Outside 10R or nonfinite target');
ending = 'END_NEWTON_RADIUS_FAILURE';
reasons = [0 1 2];
if strcmp(version,'NEWTON_PROJECTION_FAILURE_V2')
    report.failure_kind = 'newton';
    report.failure_reason = value('failure_reason');
    ending = 'END_NEWTON_PROJECTION_FAILURE';
    reasons = [0 1 2 3];
end
report.stage = value('stage');
if ~ismember(report.stage,{'vertices','rv_nodes'})
    error('edgepreserve:newtonDiagnostic','Invalid diagnostic stage.');
end
report.refinement = number('refinement');
report.iteration = number('iteration');
report.center = sscanf(value('center'),'%f');
report.radius = number('radius');
report.limit_multiplier = number('limit_multiplier');
npoints = number('num_points');
nrejected = number('num_rejected');
validateattributes(report.center,{'double'},{'numel',3,'finite'});
validateattributes(report.radius,{'double'},{'scalar','finite','positive'});
validateattributes(npoints,{'double'},{'scalar','integer','nonnegative'});
validateattributes(nrejected,{'double'},{'scalar','integer','nonnegative','<=',npoints});
expected = ['point_index base_x base_y base_z target_x target_y target_z ' ...
    'rejected reason converged step_available step_length step_iteration distance_over_radius'];
if ~strcmp(value('columns'),expected)
    error('edgepreserve:newtonDiagnostic','Unexpected diagnostic columns.');
end
body = strtrim(fread(fid,inf,'*char').');
if ~endsWith(body,ending)
    error('edgepreserve:newtonDiagnostic','Missing diagnostic end marker.');
end
body = body(1:end-numel(ending));
% gfortran writes "Infinity"; MATLAB's numeric scanner recognizes "Inf".
body = strrep(body,'Infinity','Inf');
rows = textscan(body,repmat('%f',1,14),npoints,'CollectOutput',true);
data = rows{1};
if ~isequal(size(data),[npoints 14]) || ...
        ~isequal(data(:,1),(1:npoints).') || ...
        any(~ismember(data(:,8),[0 1])) || any(~ismember(data(:,9),reasons)) || ...
        any(~ismember(data(:,10:11),[0 1]),'all') || ...
        nnz(data(:,8)) ~= nrejected || any(logical(data(:,8)) ~= (data(:,9) ~= 0))
    error('edgepreserve:newtonDiagnostic','Incomplete or inconsistent diagnostic rows.');
end
report.index = data(:,1);
report.base = data(:,2:4).';
report.target = data(:,5:7).';
report.rejected = logical(data(:,8));
report.reason = data(:,9);  % 1: outside sphere; 2: nonfinite; 3: Newton stop.
report.converged = logical(data(:,10));
report.step_available = logical(data(:,11));
report.step_length = data(:,12);
report.step_iteration = data(:,13);
report.distance_over_radius = data(:,14);

    function text = value(key)
        line = fgetl(fid);
        prefix = [key ' '];
        if ~ischar(line) || ~startsWith(line,prefix)
            error('edgepreserve:newtonDiagnostic','Expected diagnostic field %s.',key);
        end
        text = strtrim(line(numel(prefix)+1:end));
    end
    function x = number(key)
        x = str2double(value(key));
        validateattributes(x,{'double'},{'scalar','finite'});
    end
end

function [report, filename] = capture_newton_failure(exception, destination)
%CAPTURE_NEWTON_FAILURE Load and optionally preserve a controlled Newton failure.
report = [];
filename = '';
if ~ismember(exception.identifier,{'MULTISCALE_MESHER:NewtonRadiusGuard', ...
        'MULTISCALE_MESHER:NewtonProjectionFailure'}), return; end
source = diagnostic_path(exception);
try
    if isempty(source)
        error('edgepreserve:newtonDiagnostic','Newton diagnostic path is unavailable.');
    end
    report = surfsmooth3d.edgepreserve.read_newton_failure(source);
    filename = source;
    if nargin > 1 && ~isempty(destination)
        [saved,message] = copyfile(source,destination);
        if saved
            filename = destination;
            report.filename = destination;
        else
            filename = '';
            warning('edgepreserve:newtonDiagnostic','Cannot preserve Newton snapshot: %s',message);
        end
    end
catch diagnosticError
    warning('edgepreserve:newtonDiagnostic', ...
        'Newton stopped safely, but diagnostic capture failed: %s',diagnosticError.message);
end
end

function filename = diagnostic_path(exception)
filename = '';
for k = 1:numel(exception.cause)
    cause = exception.cause{k};
    if strcmp(cause.identifier,'MULTISCALE_MESHER:NewtonDiagnosticFile')
        filename = cause.message;
        return
    end
    filename = diagnostic_path(cause);
    if ~isempty(filename), return; end
end
end

function check_error(ier,root,stage)
if ier==0, return; end
id = 'adaptiveblend:nativeFailure';
message = sprintf('%s failed (Fortran ier=%d).',stage,ier);
if ismember(ier,[3 4 7])
    id = 'MULTISCALE_MESHER:NewtonProjectionFailure';
    if ier==7, id = 'MULTISCALE_MESHER:NewtonRadiusGuard'; end
    message = [message ' No incomplete child batch was accepted.'];
elseif ier==6
    message = 'Two-stage vertex projection produced a degenerate scaffold.';
elseif ier==8
    message = 'The requested mesh exceeds max_points.';
end
exception = MException(id,'%s',message);
report = [root '_newton_failure.txt'];
if isfile(report)
    preserved = [tempname '_newton_failure.txt'];
    [ok,msg] = copyfile(report,preserved);
    if ok
        exception = addCause(exception,MException('MULTISCALE_MESHER:NewtonDiagnosticFile','%s',preserved));
    else
        warning('adaptiveblend:diagnostic','Could not preserve Newton report: %s',msg);
    end
end
throw(exception);
end

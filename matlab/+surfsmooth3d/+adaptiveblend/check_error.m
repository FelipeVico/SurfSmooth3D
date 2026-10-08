function check_error(ier,root,stage)
if ier==0, return; end
id = 'adaptiveblend:nativeFailure';
message = sprintf('%s failed (Fortran ier=%d).',stage,ier);
if ismember(ier,[3 4 7])
    id = 'MULTISCALE_MESHER:NewtonProjectionFailure';
    if ier==7
        id = 'MULTISCALE_MESHER:NewtonRadiusGuard';
        if isfile([root '_newton_recovery.txt'])
            message = [message ' Deferred recovery found no acceptable local outward crossing.'];
        end
    end
    message = [message ' No incomplete child batch was accepted.'];
elseif ier==6
    message = 'Projected scaffold or recovered-patch validation produced an invalid geometry.';
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
exception = surfsmooth3d.attach_recovery_report(exception,[root '_newton_recovery.txt']);
throw(exception);
end

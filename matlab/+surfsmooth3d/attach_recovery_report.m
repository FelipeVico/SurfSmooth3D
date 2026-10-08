function exception = attach_recovery_report(exception,filename)
%ATTACH_RECOVERY_REPORT Keep a recovery report before its work folder is removed.
if ~isfile(filename), return; end
preserved = [tempname '_newton_recovery.txt'];
[ok,message] = copyfile(filename,preserved);
if ok
    exception = addCause(exception,MException( ...
        'MULTISCALE_MESHER:NewtonRecoveryFile','%s',preserved));
else
    warning('MULTISCALE_MESHER:RecoveryDiagnosticFile', ...
        'Could not preserve recovery diagnostic %s: %s',filename,message);
end
end

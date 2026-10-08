function test_newton_recovery_report()
%TEST_NEWTON_RECOVERY_REPORT Preserve stage metadata and rejected nonfinite steps.
directory = tempname; mkdir(directory);
cleanup = onCleanup(@() rmdir(directory,'s'));
filename = fullfile(directory,'recovery.txt');
empty = surfsmooth3d.read_newton_recovery(filename,false);
assert(~empty.enabled && ~empty.report_available && empty.deferred==0 && isempty(empty.events));
assert(surfsmooth3d.validate_newton_recovery(true));
assert(~surfsmooth3d.validate_newton_recovery(0));
for value = {NaN,Inf,2,-1,[true false],1+1i,'true'}
    must_fail(@() surfsmooth3d.validate_newton_recovery(value{1}));
end
fid = fopen(filename,'w');
fprintf(fid,['NEWTON_RECOVERY_V1\nBEGIN 1 0 1\nTARGET 7 2 4\n', ...
    'OUTCOME recovered\nBASE_POINT 1 2 3\nNORMAL 0 0 2\nINITIAL_POINT 1 2 3\nLAST_POINT 1 2 3.2\nINITIAL_HEIGHT 0\nLAST_HEIGHT 1D-1\n', ...
    'REJECTED_HEIGHT Infinity\nREJECTED_POINT NaN Infinity -Infinity\n', ...
    'LAST_RESIDUAL 2D-2\nLAST_GRADIENT 0 0 -1\nSIGMA 0.1\n', ...
    'INTERVAL -0.8 0.8\nROOT 0.12\nRESIDUAL 1e-10\nSLOPE -0.5\n', ...
    'POSITION_TOLERANCE 1e-10\nCOUNTS 271 1 14\n', ...
    'BRACKET 0.1 0.2 0.02 -0.08\nCANDIDATE 0.12 1e-10 -0.5 14\nCANDIDATE_OUTCOME acceptable\nEND_TARGET\nEND 1 0\n', ...
    'BEGIN 2 3 1\nTARGET 9 1 2\nOUTCOME no_bracket\n', ...
    'END_TARGET\nEND 0 1\nVALIDATION_FAILED folded_patch\n']);
fclose(fid);
info = surfsmooth3d.read_newton_recovery(filename,true,true);
assert(info.enabled && info.report_available && info.deferred==2 && info.recovered==1 && info.unresolved==1);
assert(numel(info.events)==2 && info.events(2).refinement==3 && info.validation_failed);
target = info.events(1).targets(1);
assert(target.index==7 && target.reason==2 && strcmp(target.outcome,'recovered'));
assert(isequal(target.base_point,[1 2 3]) && isequal(target.normal,[0 0 2]));
assert(abs(target.last_height-.1)<eps && isinf(target.rejected_height));
assert(isnan(target.rejected_point(1)) && target.rejected_point(3)==-Inf);
assert(isequal(target.counts,[271 1 14]) && isequal(size(target.brackets),[1 4]));
assert(isfile(info.report_file) && ~strcmp(info.report_file,filename));
preservedCleanup = onCleanup(@() delete(info.report_file));
exception = surfsmooth3d.attach_recovery_report(MException('test:failure','failed'),filename);
assert(isscalar(exception.cause) && strcmp(exception.cause{1}.identifier, ...
    'MULTISCALE_MESHER:NewtonRecoveryFile'));
exceptionReport = exception.cause{1}.message;
exceptionCleanup = onCleanup(@() delete(exceptionReport));
delete(filename);
assert(isfile(info.report_file) && isfile(exceptionReport));
fid = fopen(filename,'w'); fprintf(fid,'NEWTON_RECOVERY_V1\nBEGIN 1 0 1\n'); fclose(fid);
must_fail(@() surfsmooth3d.read_newton_recovery(filename));
fprintf('PASS: recovery switch validation, tagged reports, stage totals, nonfinite diagnostics and preservation.\n');
end

function must_fail(call)
try
    call();
catch exception
    assert(startsWith(exception.identifier,'MULTISCALE_MESHER:'));
    return
end
error('recovery_test:expectedFailure','Expected invalid input/report to fail.');
end

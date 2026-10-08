function run_surfsmooth3d_tests(mode)
%RUN_SURFSMOOTH3D_TESTS Run the extracted smoother regressions.
% mode = 'pure' checks MATLAB math without native binaries.
% mode = 'native' checks the three smoother MEX interfaces.
% mode = 'gui' checks hidden GUI workflows and independent accuracy.
% mode = 'diagnostics' reads fixtures produced by make test-native.
% mode = 'all' (default) runs all four groups. STEP has its own tests.
if nargin < 1, mode = 'all'; end
mode = validatestring(mode,{'pure','native','gui','diagnostics','all'});
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
addpath(fullfile(root,'tests','matlab'));
package = fullfile(root,'tests','fixtures','cad', ...
    'convergence_intersection_two_balls','p4_mf0p10000_rl0');
output = fullfile(root,'build','matlab-test-output');
if ~isfolder(output), mkdir(output); end
if ismember(mode,{'pure','all'})
    test_edgepreserve_fixed_cad;
    test_adaptive_reference_maps;
    test_adaptive_blend_math;
    test_newton_recovery_report;
    check_adaptive_code;
end
if ismember(mode,{'native','all'})
    test_adaptive_smoother;
    test_adaptive_blend;
    test_adaptive_blend('limits');
    test_sigma_modes(package);
    test_newton_recovery_switch;
end
if ismember(mode,{'gui','all'})
    runOutput = tempname(output); mkdir(runOutput);
    smoke_edgepreserve_gui_batch(fullfile(runOutput,'edge-sphere'));
    smoke_edgepreserve_real_package([],fullfile(runOutput,'edge-cad'));
    smoke_sigma_modes_gui(package,fullfile(runOutput,'sigma'));
    smoke_adaptive_gui(package,fullfile(runOutput,'adaptive'));
    test_adaptive_accuracy(fullfile(runOutput,'adaptive','adaptive_real_result.mat'),package);
    smoke_adaptive_blend_gui(package,fullfile(runOutput,'blend'));
    test_adaptive_blend_accuracy(fullfile(runOutput,'blend','result.mat'));
end
if ismember(mode,{'diagnostics','all'})
    report = fullfile(root,'build','native','test-output','guard');
    assert(isfile([report '_snapshot_newton_failure.txt']), ...
        'Run make test-native first to produce the native Newton diagnostic fixtures.');
    test_newton_radius_guard(report,fullfile(output,'radius'));
    test_newton_projection_failure(report,fullfile(output,'projection'));
    test_newton_guard_batch(package,[report '_snapshot_newton_failure.txt'], ...
        fullfile(output,'radius-batch'));
    test_newton_guard_batch(package,[report '_newton_newton_failure.txt'], ...
        fullfile(output,'projection-batch'));
end
fprintf('PASS: SurfSmooth3D MATLAB group %s.\n',mode);
end

function report = test_gateway_coexistence(bieRoot, outputDir)
%TEST_GATEWAY_COEXISTENCE Exercise all six gateways in one MATLAB process.
% Run in a fresh MATLAB process. Numerical routines and original packages
% are not modified; each wrapper receives its own staged scaffold directory.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
restoredefaultpath;
oldMatlab = fullfile(bieRoot,'matlab');
newMatlab = fullfile(root,'matlab');
addpath(oldMatlab,fullfile(oldMatlab,'src'),newMatlab);
if ~isfolder(outputDir), mkdir(outputDir); end
names = {'fmm3dbie_routs','fmm3dbie_adaptive_routs','fmm3dbie_adaptive_blend_routs', ...
    'surfsmooth3d_routs','surfsmooth3d_adaptive_routs','surfsmooth3d_adaptive_blend_routs'};
report = struct('matlab_version',version,'cases',struct(),'gateway_paths',struct());
for k=1:numel(names)
    expected = oldMatlab;
    if k>3, expected = newMatlab; end
    resolved = which(names{k});
    assert(startsWith(resolved,[expected filesep]),'Unexpected gateway: %s',resolved);
    report.gateway_paths.(names{k})=resolved;
end

% Load and evaluate the original BIE gateway before any smoother state exists.
S = surfer.load_from_file(fullfile(root,'tests','fixtures','spheres','sphere_unit_na2_p8.go3'));
assert(strcmp(class(S),'surfer'));
targets = struct('r',[2 -1.8 .5;.3 .2 -2.4;.7 -.6 1.1]);
density = S.r(3,:).'./vecnorm(S.r,2,1).';
exact = targets.r(3,:).'./(3*vecnorm(targets.r,2,1).'.^3);
before = lap3d.dirichlet.eval(S,density,targets,1e-8,[1 0]);
assert(norm(before-exact)/norm(exact)<2e-6);
report.bie = struct('nodes',S.npts,'patches',S.npatches, ...
    'density','z on the unit sphere','relative_error_before',norm(before-exact)/norm(exact));

for pass=1:2
    if pass==1
        addpath(oldMatlab,fullfile(oldMatlab,'src')); addpath(newMatlab);
        label='new_root_first';
    else
        addpath(newMatlab); addpath(oldMatlab,fullfile(oldMatlab,'src'));
        label='old_root_first';
    end
    assert(startsWith(which('surfer'),[oldMatlab filesep]));
    assert(startsWith(which('surfsmooth3d.surfer'),[newMatlab filesep]));
    assert(startsWith(which('multiscale_mesher'),[oldMatlab filesep]));
    assert(startsWith(which('surfsmooth3d.multiscale_mesher'),[newMatlab filesep]));
    for k=1:numel(names)
        assert(strcmp(which(names{k}),report.gateway_paths.(names{k})));
    end

    work = fullfile(outputDir,label); mkdir(work);
    [oldMesh,oldCad] = stage_inputs(root,fullfile(work,'old'));
    [newMesh,newCad] = stage_inputs(root,fullfile(work,'new'));
    opts=struct('nquad',8,'rlam',5,'adapt_sigma',1,'nrefine',0, ...
        'two_stage_smoother',false,'filetype',3,'max_refine',0, ...
        'eps_adapt',.99,'max_points',200000);
    oldOpts=opts; oldOpts.fcad=oldCad;
    newOpts=opts; newOpts.fcad=newCad;

    newA=surfsmooth3d.multiscale_mesher(newMesh,4,newOpts);
    oldA=multiscale_mesher(oldMesh,4,oldOpts);
    entry=struct('one_step',compare_surfaces(oldA{1},newA{1}));
    [oldAdaptive,oldInfo]=multiscale_mesher_adaptive(oldMesh,4,oldOpts);
    [newAdaptive,newInfo]=surfsmooth3d.multiscale_mesher_adaptive(newMesh,4,newOpts);
    entry.adaptive=compare_surfaces(oldAdaptive,newAdaptive);
    assert(isequal(oldInfo.maps,newInfo.maps));
    assert(isequal(oldInfo.parent_ids,newInfo.parent_ids));
    assert(isequal(oldInfo.depth,newInfo.depth));
    entry.adaptive.indicator_max_difference=max(abs(oldInfo.indicators-newInfo.indicators),[],'all');
    assert(entry.adaptive.indicator_max_difference<2e-11);

    % Both independently owned blend sessions are live at the same time.
    [oldHandle,ier]=fmm3dbie_adaptive_blend_routs('open',oldMesh,oldCad, ...
        fullfile(work,'old_blend'),8,4,1,5,200000); assert(ier==0);
    oldCleanup=onCleanup(@() fmm3dbie_adaptive_blend_routs('close',oldHandle));
    [oldValues,oldMaps,oldSphere]=fmm3dbie_adaptive_blend_routs('get',oldHandle);
    [newHandle,ier]=surfsmooth3d_adaptive_blend_routs('open',newMesh,newCad, ...
        fullfile(work,'new_blend'),8,4,1,5,200000); assert(ier==0);
    newCleanup=onCleanup(@() surfsmooth3d_adaptive_blend_routs('close',newHandle));
    [newValues,newMaps,newSphere]=surfsmooth3d_adaptive_blend_routs('get',newHandle);
    [oldAgain,oldMapsAgain,oldSphereAgain]=fmm3dbie_adaptive_blend_routs('get',oldHandle);
    assert(isequal(oldValues,oldAgain) && isequal(oldMaps,oldMapsAgain) && isequal(oldSphere,oldSphereAgain));
    assert(isequal(oldMaps,newMaps) && isequal(oldSphere,newSphere));
    entry.blend=struct('both_sessions_live',true,'old_session_unchanged',true, ...
        'max_difference',max(abs(oldValues-newValues),[],'all'));
    assert(entry.blend.max_difference<2e-11);
    entry.loaded_mex=check_loaded(names,report.gateway_paths);
    clear oldCleanup;
    [afterClose,afterMaps,afterSphere]=surfsmooth3d_adaptive_blend_routs('get',newHandle);
    assert(isequal(newValues,afterClose) && isequal(newMaps,afterMaps) && isequal(newSphere,afterSphere));
    entry.blend.new_session_survives_old_close=true;
    clear newCleanup;

    % A changed original-library call must not contaminate the new library.
    oldOpts.rlam=4;
    oldB=multiscale_mesher(oldMesh,4,oldOpts);
    newAgain=surfsmooth3d.multiscale_mesher(newMesh,4,newOpts);
    entry.interleaved_A_B_A=compare_surfaces(newA{1},newAgain{1});
    a=extract_arrays(newA{1}); b=extract_arrays(oldB{1});
    entry.changed_parameter_max_difference=max(abs(a-b),[],'all');
    assert(entry.changed_parameter_max_difference>1e-8);
    after=lap3d.dirichlet.eval(S,density,targets,1e-8,[1 0]);
    entry.bie_max_difference_from_before=max(abs(after-before));
    entry.bie_relative_error=norm(after-exact)/norm(exact);
    assert(entry.bie_max_difference_from_before<2e-12 && entry.bie_relative_error<2e-6);
    report.cases.(label)=entry;
    fprintf('GATEWAY_COEXISTENCE_PATH_ORDER_PASSED: %s\n',label);
end
report.all_six_loaded=check_loaded(names,report.gateway_paths);
clear mex;
[~,remaining]=inmem('-completenames');
report.mex_after_clear=remaining;
report.clear_mex_returned=true;
report.normal_exit='Checked by the invoking fresh MATLAB process exit status.';
fid=fopen(fullfile(outputDir,'coexistence-report.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('GATEWAY_COEXISTENCE_PASSED: all six gateways, two path orders, live sessions, A/B/A and original BIE evaluation.\n');
end

function [mesh,cad]=stage_inputs(root,directory)
mkdir(directory);
mesh=fullfile(directory,'cylinder.gidmsh');
cad=fullfile(directory,'cylinder_skeleton.txt');
copyfile(fullfile(root,'tests','fixtures','cylinder','cylinder.gidmsh'),mesh);
copyfile(fullfile(root,'tests','fixtures','cylinder','cylinder_skeleton.txt'),cad);
end

function report=compare_surfaces(a,b)
[av,ac,ao,ai,at,aw]=extract_arrays(a);
[bv,bc,bo,bi,bt,bw]=extract_arrays(b);
assert(isequal(ao,bo) && isequal(ai,bi) && isequal(at,bt));
report=struct('patches',a.npatches,'nodes',a.npts, ...
    'values_max_difference',max(abs(av-bv),[],'all'), ...
    'coefficients_max_difference',max(abs(ac-bc),[],'all'), ...
    'weights_max_difference',max(abs(aw-bw),[],'all'), ...
    'values_exact',isequal(av,bv),'coefficients_exact',isequal(ac,bc),'weights_exact',isequal(aw,bw));
assert(report.values_max_difference<2e-11 && report.coefficients_max_difference<2e-11 && report.weights_max_difference<2e-11);
end

function paths=check_loaded(names,expected)
[~,paths]=inmem('-completenames');
for k=1:numel(names)
    assert(any(strcmp(paths,expected.(names{k}))),'Gateway not in memory: %s',names{k});
end
paths=paths(ismember(paths,struct2cell(expected)));
assert(numel(paths)==6);
end

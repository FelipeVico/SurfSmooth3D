function capture_equivalence(mode, legacyRoot, packageDir, outputDir)
%CAPTURE_EQUIVALENCE Capture numerical results without serialized class objects.
% Run each mode in a fresh MATLAB process, then use compare_equivalence.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
restoredefaultpath;
if strcmp(mode,'reference')
    addpath(fullfile(legacyRoot,'matlab'),fullfile(legacyRoot,'matlab','src'));
    prefix = '';
elseif strcmp(mode,'extracted')
    addpath(fullfile(root,'matlab'));
    setup_surfsmooth3d;
    prefix = 'surfsmooth3d.';
    assert(isempty(which('fmm3dbie_routs')),'BIE gateway on standalone path.');
else
    error('Use reference or extracted mode.');
end
if ~isfolder(outputDir), mkdir(outputDir); end
api = @(name) str2func([prefix name]);
loadPackage = api('edgepreserve.load_package');
stage = api('edgepreserve.stage_source');
settingsFcn = api('edgepreserve.default_settings');
mesherOpts = api('edgepreserve.mesher_options');
solve = api('multiscale_mesher');
sigmaEval = api('multiscale_mesher_sigma_eval');
writeGo3 = api('edgepreserve.write_go3');
pkg = loadPackage(packageDir);
[work,cleanup] = stage(pkg); %#ok<ASGLU>
settings = settingsFcn();
settings.selectedEdgeIds = pkg.edges(1).edge_id;
opts = mesherOpts(pkg,work,settings,1);
opts.two_stage_smoother = false;
one = solve(work.scaffold_file,4,opts);
record = struct();
record.source = surface_arrays(pkg.surface);
record.edges = pkg.edges;
record.one_step = cellfun(@surface_arrays,one,'UniformOutput',false);
writeGo3(fullfile(outputDir,'one_step.go3'),one{end});
opts.two_stage_smoother = true;
two = solve(work.scaffold_file,4,opts);
record.two_stage = cellfun(@surface_arrays,two,'UniformOutput',false);
writeGo3(fullfile(outputDir,'two_stage.go3'),two{end});
targets = pkg.surface.r(:,1:37:end);
for modeIndex = 0:3
    sigmaOpts = opts; sigmaOpts.adapt_sigma = modeIndex;
    [sigma,gradient] = sigmaEval(work.scaffold_file,4,sigmaOpts,targets);
    record.sigma{modeIndex+1} = struct('value',sigma,'gradient',gradient);
end
pairFcn = api('edgepreserve.compute_pair');
pair = pairFcn(pkg,work,two{2},1,settings,sigmaEval);
assert(isempty(pair.blend_error),pair.blend_error);
assert(pair.report_smooth.is_valid && pair.report_blend.is_valid);
record.uniform_blend = struct('surface',surface_arrays(pair.S_blend), ...
    'reference',surface_arrays(pair.S_cad),'beta',pair.beta,'dsoft',pair.dsoft, ...
    'sigma',pair.sigma,'gradient',pair.gradSigma);
writeGo3(fullfile(outputDir,'uniform_blend.go3'),pair.S_blend);
adaptiveOpts = struct('order',4,'rlam',2,'adapt_sigma',1, ...
    'eps_adapt',.01,'max_refine',1,'max_points',200000);
adaptive = api('adaptivesmoother.solve');
a = adaptive(pkg,adaptiveOpts);
record.adaptive = struct('surface',surface_arrays(a.surface), ...
    'reference',surface_arrays(a.reference),'info',numeric_info(a.info));
writeGo3(fullfile(outputDir,'adaptive.go3'),a.surface);
blendSettings = api('adaptiveblend.default_settings');
s = blendSettings(); s.order = 4; s.max_refine = 1; s.eps_adapt = .01;
s.selectedEdgeIds = settings.selectedEdgeIds;
blend = api('adaptiveblend.solve');
b = blend(pkg,s);
record.adaptive_blend = struct('surface',surface_arrays(b.surface), ...
    'smooth',surface_arrays(b.smooth),'reference',surface_arrays(b.reference), ...
    'info',numeric_info(b.info),'beta',b.beta,'sigma',b.sigma);
writeGo3(fullfile(outputDir,'adaptive_blend.go3'),b.surface);
% Repeat the initial calculation after all saved-state consumers have run.
opts.two_stage_smoother = false;
again = solve(work.scaffold_file,4,opts);
record.repeat = cellfun(@surface_arrays,again,'UniformOutput',false);
assert(isequaln(record.one_step,record.repeat),'Repeated A-B-A calculation changed.');
save(fullfile(outputDir,'numerical-results.mat'),'record','-v7');
fprintf('CAPTURE_EQUIVALENCE_PASSED: %s\n',mode);
clear mex;
end

function result = surface_arrays(S)
[values,coefs,orders,offsets,types,weights] = extract_arrays(S);
result = struct('values',values,'coefs',coefs,'orders',orders, ...
    'offsets',offsets,'types',types,'weights',weights);
end

function result = numeric_info(info)
names = {'parent_ids','depth','maps','indicators','independent_checked', ...
    'unresolved','unresolved_count','converged','stop_reason','achieved_depth', ...
    'launch_diameter','center','radius'};
result = struct();
for k = 1:numel(names)
    if isfield(info,names{k}), result.(names{k}) = info.(names{k}); end
end
result.history = table2array(info.history);
end

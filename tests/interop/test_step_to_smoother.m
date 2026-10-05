function report = test_step_to_smoother(installationRoot, sourceRoot, outputDir)
%TEST_STEP_TO_SMOOTHER Exercise a fresh STEP package through all smoothers.
% sourceRoot supplies test inputs only; all implementation resolves in install.
restoredefaultpath;
addpath(fullfile(installationRoot,'matlab'));
setup_surfsmooth3d;
assert(isempty(which('fmm3dbie_routs')));
assert(isempty(which('surfer')));
assert(startsWith(which('surfsmooth3d.multiscale_mesher'),installationRoot));
assert(startsWith(which('surfsmooth3d.stepmesher.mesh_step'),installationRoot));
if ~isfolder(outputDir), mkdir(outputDir); end
filename=fullfile(sourceRoot,'step_mesher','examples','step_files','intersection_two_balls.step');
opts=struct('order',4,'edge_gll_order',8,'mesh_fraction',.1, ...
    'refinement_level',0,'mesh_profile','rigid');
mesh=surfsmooth3d.stepmesher.mesh_step(filename,opts);
assert(strcmp(mesh.meta.cad_identity,'native_occ_topology_v1'));
packageDir=fullfile(outputDir,'cad');
surfsmooth3d.stepmesher.export_geom_package(mesh,packageDir,'intersection_two_balls');
pkg=surfsmooth3d.edgepreserve.load_package(packageDir);
settings=surfsmooth3d.edgepreserve.default_settings();
settings.orders=4; settings.levels=0;
settings.selectedEdgeIds=pkg.edges(1).edge_id;
[manifest,runDir]=surfsmooth3d.edgepreserve.run_batch(pkg,settings,fullfile(outputDir,'uniform'));
assert(all(manifest.status=="saved"));
aopts=struct('order',4,'rlam',2,'adapt_sigma',1,'eps_adapt',.01, ...
    'max_refine',1,'max_points',200000);
result=surfsmooth3d.adaptivesmoother.solve(pkg,aopts);
assert(result.report.is_valid);
adaptiveFile=surfsmooth3d.adaptivesmoother.save_result(pkg,result,fullfile(outputDir,'adaptive'));
bopts=surfsmooth3d.adaptiveblend.default_settings();
bopts.order=4; bopts.max_refine=1; bopts.eps_adapt=.01;
bopts.selectedEdgeIds=pkg.edges(1).edge_id;
blend=surfsmooth3d.adaptiveblend.solve(pkg,bopts);
assert(blend.report.is_valid && blend.report_smooth.is_valid);
blendManifest=surfsmooth3d.adaptiveblend.save_result(pkg,blend,fullfile(outputDir,'blend'));
assert(all(blendManifest.status=="saved"));
for k=1:numel(pkg.source_paths)
    assert(strcmp(pkg.source_hashes{k},surfsmooth3d.edgepreserve.file_sha256(pkg.source_paths{k})));
end
report=struct('installed_root',installationRoot,'step_patches',mesh.npatches, ...
    'uniform_exports',height(manifest),'uniform_directory',runDir, ...
    'fully_smooth_file',fullfile(runDir,manifest.filename(1)), ...
    'adaptive_patches',result.surface.npatches,'adaptive_converged',result.info.converged, ...
    'adaptive_file',adaptiveFile,'blend_patches',blend.surface.npatches, ...
    'blend_converged',blend.info.converged,'blend_exports',height(blendManifest), ...
    'cad_sources_unchanged',true,'bie_on_path',false);
fid=fopen(fullfile(outputDir,'pipeline-report.json'),'w');
c=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('STEP_TO_SMOOTHER_PASSED: %d STEP patches, all four workflow exports validated.\n',mesh.npatches);
clear mex;
end

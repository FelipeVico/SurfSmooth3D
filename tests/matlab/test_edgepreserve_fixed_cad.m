function test_edgepreserve_fixed_cad()
%TEST_EDGEPRESERVE_FIXED_CAD Reference maps, derivatives, blending and exports.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
for sourceLevel = 0:2
    pkg = polynomial_package(sourceLevel);
    for outputLevel = 0:3
        for order = [2 4 6]
            output = polynomial_surface(outputLevel,order,'smoother',false);
            [cad,map] = surfsmooth3d.edgepreserve.resample_cad(pkg,output,outputLevel);
            expected = polynomial_positions(map.coarse_uv,2);
            assert(max(abs(cad.r(:)-expected(:))) < 1e-11*max(abs(expected(:))));
        end
    end
end
% Boundary ties must go to the first matching source triangle.
vertices = [0 1 0 .5 .5 0;0 0 1 0 .5 .5];
map = surfsmooth3d.edgepreserve.reference_map(1,0,1,vertices);
assert(isequal(map.source_patch,[1 4 3 1 2 1]));
children = surfsmooth3d.edgepreserve.triangle_maps(1,'smoother');
assert(isequal(children(:,:,2),[.5 -.5 0;.5 0 -.5]));
assert(all(arrayfun(@(k) det(children(:,2:3,k))>0,1:4)));
didFail = false;
try
    surfsmooth3d.edgepreserve.reference_map(1,0,1,[2;0]);
catch
    didFail = true;
end
assert(didFail);
for level = 0:3
    planar = polynomial_surface(level,4,'smoother',true);
    jac = vecnorm(cross(planar.du,planar.dv,1));
    assert(max(abs(jac-4^-level)) < 1e-12);
    pkg = polynomial_package(0,true);
    matching = surfsmooth3d.edgepreserve.resample_cad(pkg,planar,level);
    report = surfsmooth3d.edgepreserve.validity_report(planar,'plane',matching);
    assert(report.is_valid && max(abs(report.jac_ratio_to_ref-1)) < 1e-11);
end
assert(isequal(surfsmooth3d.edgepreserve.parse_integer_list('[8,4 6 4]',1,20),[4 6 8]));
assert(isempty(surfsmooth3d.edgepreserve.parse_integer_list('[]',1,inf,true)));
for invalid = {'4;system','[1:3]','1.5','-1','21'}
    didFail = false;
    try
        surfsmooth3d.edgepreserve.parse_integer_list(invalid{1},1,20);
    catch
        didFail = true;
    end
    assert(didFail);
end
test_blend();
test_batch_failures();
fprintf('PASS: fixed-CAD mapping, central child, boundary ties, derivatives, blend, batch/export failures.\n');
end

function pkg = polynomial_package(level,planar)
if nargin < 2, planar = false; end
pkg = struct('surface',polynomial_surface(level,4,'step',planar), ...
    'source_order',4,'source_level',level,'coarse_triangle_count',2);
end

function S = polynomial_surface(level,order,convention,planar)
maps = surfsmooth3d.edgepreserve.triangle_maps(level,convention);
uv = surfsmooth3d.internal.koorn.rv_nodes(order);
coarse = zeros(2,size(uv,2)*size(maps,3));
for k = 1:size(maps,3)
    coarse(:,(k-1)*size(uv,2)+(1:size(uv,2))) = maps(:,1,k)+maps(:,2:3,k)*uv;
end
r = polynomial_positions(coarse,2,planar);
npatches = 2*4^level;
n = size(uv,2);
template = struct('npatches',npatches,'npts',n*npatches, ...
    'norders',order*ones(npatches,1),'iptype',ones(npatches,1), ...
    'r',r,'ixyzs',(1:n:n*npatches+1).');
src = surfsmooth3d.edgepreserve.srcvals_from_positions(template,r);
S = surfsmooth3d.surfer(npatches,template.norders,src,template.iptype);
end

function r = polynomial_positions(uv,parents,planar)
if nargin < 3, planar = false; end
r = [];
for parent = 1:parents
    u = uv(1,:); v = uv(2,:);
    z = parent+.1*u.^2+.07*u.*v+.08*v.^2;
    if planar, z = parent*ones(size(u)); end
    r = [r [u+10*parent;v+parent;z]]; %#ok<AGROW>
end
end

function test_blend()
pkg = polynomial_package(0);
cad = pkg.surface;
smooth = surfsmooth3d.surfer(cad.npatches,cad.norders, ...
    surfsmooth3d.edgepreserve.srcvals_from_positions(cad,cad.r+[0;0;.15]),cad.iptype);
edgeX = [linspace(9,12,50);ones(1,50);ones(1,50)];
panel = struct('xyz',edgeX,'w',ones(1,50)*3/49);
edges = struct('edge_id',1,'panels',panel);
params = surfsmooth3d.edgepreserve.default_settings();
params.sigma = ones(cad.npts,1);
params.geometryDiameter = 15;
params.selectedEdgeIds = 1;
[blend,beta,dsoft] = surfsmooth3d.edgepreserve.blend(cad,smooth,edges,params);
dist2 = sum((reshape(cad.r,3,[],1)-reshape(edgeX,3,1,[])).^2,1);
A = reshape(exp(-dist2/(2*params.c_eps^2)),cad.npts,[])*panel.w(:) / ...
    (sqrt(2*pi)*params.c_eps);
expectedD = params.c_eps*sqrt(max(0,-2*log(max(A,params.distanceFloor))));
expectedBeta = .5*erfc((expectedD-params.c_rad)/(sqrt(2)*params.c_tau));
assert(max(abs(beta-expectedBeta))<1e-13 && max(abs(dsoft-expectedD))<1e-13);
assert(min(beta)<1e-5 && max(beta)>.9);
expectedR = expectedBeta.'.*smooth.r+(1-expectedBeta.').*cad.r;
expectedSrc = surfsmooth3d.edgepreserve.srcvals_from_positions(cad,expectedR);
[actual,~,~,~,~,~] = extract_arrays(blend);
assert(max(abs(actual(:)-expectedSrc(:)))<1e-12);
params.selectedEdgeIds = [];
[emptyBlend,emptyBeta] = surfsmooth3d.edgepreserve.blend(cad,smooth,edges,params);
assert(all(emptyBeta==0) && max(abs(emptyBlend.r(:)-cad.r(:)))<1e-12);
assert(isequal(smooth.r,cad.r+[0;0;.15]));
report = surfsmooth3d.edgepreserve.validity_report(blend,'blend',cad);
assert(report.is_valid);
end

function test_batch_failures()
directory = tempname; mkdir(directory);
cleanup = onCleanup(@() rmdir(directory,'s'));
pkg = polynomial_package(0);
pkg.geometry_name = 'polynomial_fixture';
pkg.package_dir = directory;
pkg.scaffold_file = fullfile(directory,'scaffold.gidmsh');
write_fixture_scaffold(pkg.scaffold_file);
pkg.source_paths = {pkg.scaffold_file};
pkg.source_hashes = {surfsmooth3d.edgepreserve.file_sha256(pkg.scaffold_file)};
pkg.edges = struct('edge_id',1,'panels',struct('xyz',[10 11;1 1;1 1],'w',[.5 .5]));
settings = surfsmooth3d.edgepreserve.default_settings();
settings.orders = [4 6]; settings.levels = [0 1 2];
settings.selectedEdgeIds = 999; % full succeeds independently of partial.
hooks = struct('solve',@failure_solver,'sigma',@constant_sigma);
[manifest,runDir] = surfsmooth3d.edgepreserve.run_batch(pkg,settings,directory,@(~) [],hooks);
assert(height(manifest)==12);
assert(nnz(manifest.status=="saved")==4);
assert(all(manifest.status(manifest.refinement==2)=="failed"));
assert(all(manifest.status(manifest.mode=="edgepreserve")=="failed"));
assert(isfile(fullfile(runDir,'parameters.mat')) && isfile(fullfile(runDir,'batch_manifest.csv')));
assert(strcmp(surfsmooth3d.edgepreserve.file_sha256(pkg.scaffold_file),pkg.source_hashes{1}));
% A failed sigma evaluation still permits an independently valid full export.
hooks.sigma = @failure_sigma;
settings.orders = 4; settings.levels = 0; settings.selectedEdgeIds = 1;
[sigmaManifest,~] = surfsmooth3d.edgepreserve.run_batch(pkg,settings,directory,@(~) [],hooks);
assert(sigmaManifest.status(1)=="saved" && sigmaManifest.status(2)=="failed");
hooks.sigma = @constant_sigma;
bad = pkg.surface;
[src,~,~,~,~,~] = extract_arrays(bad); src(10:12,:) = -src(10:12,:);
bad = surfsmooth3d.surfer(bad.npatches,bad.norders,src,bad.iptype);
report = surfsmooth3d.edgepreserve.validity_report(bad,'bad',pkg.surface);
assert(~report.is_valid);
didFail = false;
try
    surfsmooth3d.edgepreserve.export_surface(fullfile(directory,'bad.go3'),bad,report);
catch
    didFail = true;
end
assert(didFail && ~isfile(fullfile(directory,'bad.go3')));
settings.orders = [4 6 8]; settings.levels = [1 2 3]; settings.selectedEdgeIds = [];
hooks.solve = @successful_solver;
[manifest,runDir] = surfsmooth3d.edgepreserve.run_batch(pkg,settings,directory,@(~) [],hooks);
assert(height(manifest)==18 && all(manifest.status=="saved"));
assert(numel(dir(fullfile(runDir,'*.go3')))==18);
end

function surfaces = failure_solver(~,order,opts)
if opts.nrefine>=2, error('fixture:solverFailure','Deliberate maximum-level failure.'); end
surfaces = successful_solver('',order,opts);
end

function surfaces = successful_solver(~,order,opts)
surfaces = cell(1,opts.nrefine+1);
for level = 0:opts.nrefine
    surfaces{level+1} = polynomial_surface(level,order,'smoother',false);
end
end

function [sigma,gradient] = constant_sigma(~,~,~,targets)
sigma = ones(size(targets,2),1); gradient = zeros(size(targets));
end

function [sigma,gradient] = failure_sigma(varargin) %#ok<STOUT>
error('fixture:sigmaFailure','Deliberate sigma evaluation failure.');
end

function write_fixture_scaffold(filename)
fid = fopen(filename,'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'MESH dimension 3 ElemType Triangle Nnode 3\nCoordinates\n');
fprintf(fid,'1 10 1 1\n2 11 1 1\n3 10 2 1\n4 20 2 2\n5 21 2 2\n6 20 3 2\n');
fprintf(fid,'End Coordinates\nElements\n1 1 2 3\n2 4 5 6\nEnd Elements\n');
end

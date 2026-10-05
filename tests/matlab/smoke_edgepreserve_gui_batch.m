function smoke_edgepreserve_gui_batch(outputParent)
%SMOKE_EDGEPRESERVE_GUI_BATCH Exercise the real GUI batch button: nine pairs.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
clear dir;
if nargin < 1, outputParent = tempname; end
mkdir(outputParent);
packageDir = fullfile(outputParent,'sphere_package');
mkdir(packageDir);
S = fixture_sphere(1,2,[0;0;0],10,1);
surfsmooth3d.edgepreserve.write_go3(fullfile(packageDir,'surface.go3'),S);
write_scaffold(fullfile(packageDir,'scaffold.gidmsh'),S);
write_edges(fullfile(packageDir,'edges_gll.txt'));
write_metadata(fullfile(packageDir,'metadata.txt'));
fig = surfsmooth3d.edgepreserve.launch_gui(root,packageDir,'off',outputParent);
cleanup = onCleanup(@() close_if_open(fig));
rlam = findobj(fig,'Tag','rlam'); rlam.String = '4';
rlam.Callback(rlam,[]);
selection = findobj(fig,'Tag','edges'); selection.String = '1';
selection.Callback(selection,[]);
batch = findobj(fig,'Tag','batch');
batch.Callback(batch,[]);
folders = dir(fullfile(outputParent,'sphere_fixture','batch_*'));
assert(isscalar(folders));
runDir = fullfile(folders(1).folder,folders(1).name);
manifest = readtable(fullfile(runDir,'batch_manifest.csv'), ...
    'Delimiter',',','TextType','string','VariableNamingRule','preserve');
assert(height(manifest)==18 && all(manifest.status=="saved"), ...
    strjoin(string(manifest.error_message),newline));
assert(numel(dir(fullfile(runDir,'*.go3')))==18);
for k = 1:height(manifest)
    saved = surfsmooth3d.surfer.load_from_file(fullfile(runDir,manifest.filename(k)));
    assert(saved.npatches==S.npatches*4^manifest.refinement(k));
    assert(all(saved.norders==manifest.order(k)));
end
assert(isempty(getappdata(fig,'edgepreserve_preview')));
exportapp(fig,fullfile(outputParent,'gui_batch.png'));
close(fig);
fprintf('PASS: real GUI batch [4 6 8] x [1 2 3], eighteen valid reloadable surfaces.\n%s\n',runDir);
end

function write_scaffold(filename,S)
vertices = [];
% surfsmooth3d.surfer.end_pt_verts stores (0,0),(0,1),(1,0), not the canonical vertex order.
% Snap extrapolated vertices to the na=2 sphere's known cube-grid coordinates.
coordinates = [-1 -1/sqrt(2) -1/sqrt(3) 0 1/sqrt(3) 1/sqrt(2) 1];
for k = 1:S.npatches
    approximate = S.end_pt_verts{k}(:,[1 3 2]);
    [~,index] = min(abs(approximate(:)-coordinates),[],2);
    vertices = [vertices reshape(coordinates(index),3,3)]; %#ok<AGROW>
end
[nodes,~,indices] = uniquetol(vertices.',1e-10,'ByRows',true);
tris = reshape(indices,3,[]);
fid = fopen(filename,'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'MESH dimension 3 ElemType Triangle Nnode 3\nCoordinates\n');
fprintf(fid,'%d %.16e %.16e %.16e\n',[(1:size(nodes,1)).' nodes].');
fprintf(fid,'End Coordinates\nElements\n');
fprintf(fid,'%d %d %d %d\n',[(1:size(tris,2));tris]);
fprintf(fid,'End Elements\n');
end

function write_edges(filename)
fid = fopen(filename,'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'STEP_MESHER_EDGES_GLL_V2\nedge_order 8\nnum_cad_edges 1\n');
fprintf(fid,'CAD_EDGE\nedge_id 1\ncurve_tag 1\nadjacent_faces 1 2\nnum_panels 8\n');
[a,w] = gll_order8();
for panel = 1:8
    theta = (panel-1+a)*2*pi/8;
    xyz = [cos(theta);sin(theta);zeros(size(theta))];
    tangent = [-sin(theta);cos(theta);zeros(size(theta))];
    fprintf(fid,'PANEL\npanel_id %d\nnode_tags %d %d\nnum_nodes 9\n',panel,panel,mod(panel,8)+1);
    fprintf(fid,'columns a x y z tx ty tz w\n');
    fprintf(fid,'%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n', ...
        [a;xyz;tangent;w*2*pi/8]);
    fprintf(fid,'END_PANEL\n');
end
fprintf(fid,'END_CAD_EDGE\nEND_STEP_MESHER_EDGES_GLL\n');
end

function [a,w] = gll_order8()
% Fixed nine-node Legendre GLL fixture, scaled from [-1,1] to [0,1].
x = [-1 -.8997579954114602 -.6771862795107378 -.3631174638261782 ...
    0 .3631174638261782 .6771862795107378 .8997579954114602 1];
w = [.02777777777777778 .1654953615608055 .2745387125001617 ...
    .3464285109730463 .3715192743764172 .3464285109730463 ...
    .2745387125001617 .1654953615608055 .02777777777777778]/2;
a = (x+1)/2;
end

function write_metadata(filename)
fid = fopen(filename,'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'STEP_MESHER_GEOM_PACKAGE_V1\nname sphere_package\n');
fprintf(fid,'source_step sphere_fixture.step\nsurface_file surface.go3\n');
fprintf(fid,'scaffold_file scaffold.gidmsh\nedges_file edges_gll.txt\n');
fprintf(fid,'surface_order 10\nedge_order 8\nrefinement_level 0\n');
fprintf(fid,'END_STEP_MESHER_GEOM_PACKAGE\n');
end

function close_if_open(fig)
if isgraphics(fig), close(fig); end
end

function S = fixture_sphere(radius, divisions, center, order, patchType)
%FIXTURE_SPHERE Frozen original fmm3dbie sphere geometry for standalone tests.
% These three fixtures were produced by the approved baseline; no BIE MEX is
% needed when running the extracted package's tests.
if patchType==1 && radius==1 && divisions==2 && isequal(center,[0;0;0]) && ismember(order,[8 10])
    name = sprintf('sphere_unit_na2_p%d.go3',order);
elseif patchType==1 && radius==4 && divisions==1 && isequal(center,[30;-20;10]) && order==4
    name = 'sphere_offset_na1_p4.go3';
else
    error('surfsmooth3d:testFixture','Requested sphere is not a frozen test fixture.');
end
S = surfsmooth3d.surfer.load_from_file(fullfile(surfsmooth3d.root(), ...
    'tests','fixtures','spheres',name));
end

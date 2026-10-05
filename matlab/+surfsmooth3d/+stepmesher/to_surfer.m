function S = to_surfer(mesh)
%TO_SURFER Construct an fmm3dbie surfer object from a stepmesher mesh.

if isempty(which('surfsmooth3d.surfer'))
    error('stepmesher:missingSurfer', ...
        'surfer was not found. Add SurfSmooth3D/matlab to the MATLAB path.');
end

[srcvals, norders, iptype] = surfsmooth3d.stepmesher.to_srcvals(mesh);
S = surfsmooth3d.surfer(double(mesh.npatches), norders, srcvals, iptype);
end

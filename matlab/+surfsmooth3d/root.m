function directory = root()
%ROOT Return the SurfSmooth3D source or installation root.
directory = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

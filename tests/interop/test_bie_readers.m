function report = test_bie_readers(bieRoot, outputDir)
%TEST_BIE_READERS Compare two MATLAB classes and the unchanged native BIE reader.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
restoredefaultpath;
addpath(fullfile(bieRoot,'matlab'),fullfile(bieRoot,'matlab','src'));
addpath(fullfile(root,'matlab'));
setup_surfsmooth3d;
files = dir(fullfile(outputDir,'*.go3'));
report = struct('files',struct());
for k=1:numel(files)
    filename=fullfile(outputDir,files(k).name);
    A=surfsmooth3d.surfer.load_from_file(filename);
    B=surfer.load_from_file(filename);
    [av,ac,ao,ai,at,aw]=extract_arrays(A);
    [bv,bc,bo,bi,bt,bw]=extract_arrays(B);
    assert(isequal(av,bv) && isequal(ac,bc) && isequal(aw,bw));
    assert(isequal(ao,bo) && isequal(ai,bi) && isequal(at,bt));
    fid=fopen([filename '.bie.bin'],'r'); assert(fid>=0);
    c=onCleanup(@() fclose(fid));
    dims=fread(fid,2,'int64'); assert(isequal(dims,[A.npatches;A.npts]));
    nv=fread(fid,[12 A.npts],'double');
    nc=fread(fid,[9 A.npts],'double');
    nw=fread(fid,A.npts,'double');
    assert(isequal(nv,av),'Native BIE reader changed source values.');
    coefError=max(abs(nc-ac),[],'all'); weightError=max(abs(nw(:)-aw(:)));
    assert(coefError<2e-11*max(1,max(abs(ac),[],'all')));
    assert(weightError<2e-13*max(1,max(abs(aw))));
    entry=struct('patches',A.npatches,'nodes',A.npts, ...
        'matlab_exact',true,'native_values_exact',true, ...
        'native_coefficient_max_error',coefError,'native_weight_max_error',weightError);
    report.files.(matlab.lang.makeValidName(files(k).name))=entry;
    clear c;
end
% Check both MATLAB path orders without introducing global class aliases.
for pass=1:2
    if pass==1, addpath(fullfile(bieRoot,'matlab')); else, addpath(fullfile(root,'matlab')); end
    assert(startsWith(which('surfer'),bieRoot));
    assert(startsWith(which('surfsmooth3d.surfer'),root));
    assert(startsWith(which('koorn.rv_nodes'),bieRoot));
    assert(startsWith(which('surfsmooth3d.internal.koorn.rv_nodes'),root));
    assert(startsWith(which('fmm3dbie_routs'),bieRoot));
    assert(startsWith(which('surfsmooth3d_routs'),root));
    A=surfsmooth3d.surfer.load_from_file(fullfile(outputDir,files(1).name));
    B=surfer.load_from_file(fullfile(outputDir,files(1).name));
    assert(isequal(extract_arrays(A),extract_arrays(B)));
end
% Supply a small closed geometry written by the independent package for the
% separate BIE-only solve. Its analytic fixture was frozen before extraction.
S=surfsmooth3d.surfer.load_from_file(fullfile(root,'tests','fixtures','spheres','sphere_unit_na2_p8.go3'));
surfsmooth3d.edgepreserve.write_go3(fullfile(outputDir,'bie_solve_input.go3'),S);
fid=fopen(fullfile(outputDir,'bie-readers-report.json'),'w');
c=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('BIE_READERS_PASSED: %d output families and both MATLAB path orders.\n',numel(files));
end

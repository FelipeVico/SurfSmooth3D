function report = compare_equivalence(referenceDir, extractedDir)
%COMPARE_EQUIVALENCE Require exact equality of the frozen same-runtime results.
a = load(fullfile(referenceDir,'numerical-results.mat'));
b = load(fullfile(extractedDir,'numerical-results.mat'));
names = fieldnames(a.record);
report = struct('exact',true,'cases',struct());
for k=1:numel(names)
    name=names{k}; same=isequaln(a.record.(name),b.record.(name));
    report.cases.(name)=same;
    assert(same,'Extraction changed numerical results for %s.',name);
end
files = dir(fullfile(referenceDir,'*.go3'));
for k=1:numel(files)
    name=files(k).name;
    a=fileread(fullfile(referenceDir,name)); b=fileread(fullfile(extractedDir,name));
    assert(strcmp(a,b),'Export differs for %s.',name);
end
report.go3_files = {files.name};
fid=fopen(fullfile(extractedDir,'equivalence-report.json'),'w');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('EXTRACTION_EQUIVALENCE_PASSED: %d exact numerical cases, %d byte-identical exports.\n',numel(names),numel(files));
end

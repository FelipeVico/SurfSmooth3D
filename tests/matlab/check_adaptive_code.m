function check_adaptive_code()
%CHECK_ADAPTIVE_CODE Check the new MATLAB path without changing old drivers.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
files = [dir(fullfile(root,'matlab','+surfsmooth3d','+adaptivesmoother','*.m')); ...
    dir(fullfile(root,'matlab','+surfsmooth3d','multiscale_mesher_adaptive.m')); ...
    dir(fullfile(root,'examples','matlab','driver_two_stage_adaptive_smooth.m')); ...
    dir(fullfile(root,'tests','matlab','*adaptive*.m'))];
count = 0;
for k = 1:numel(files)
    filename = fullfile(files(k).folder,files(k).name);
    issues = checkcode(filename,'-id');
    for j = 1:numel(issues)
        fprintf('%s:%d %s %s\n',files(k).name,issues(j).line,issues(j).id,issues(j).message);
    end
    count = count+numel(issues);
end
assert(count==0,'%d MATLAB checkcode issues.',count);
fprintf('PASS: MATLAB checkcode on adaptive additions.\n');
end

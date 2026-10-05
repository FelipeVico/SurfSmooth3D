function directory = default_input_root(repoRoot)
%DEFAULT_INPUT_ROOT Locate this package's STEP output or bundled CAD fixtures.
if nargin < 1 || isempty(repoRoot)
    repoRoot = surfsmooth3d.root();
end
candidates = {fullfile(repoRoot,'outputs','step'), ...
    fullfile(repoRoot,'tests','fixtures','cad'), char(repoRoot)};
for k = 1:numel(candidates)
    if isfolder(candidates{k})
        directory = candidates{k};
        return
    end
end
directory = char(repoRoot);
end

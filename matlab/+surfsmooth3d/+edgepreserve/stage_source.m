function [work, cleanup] = stage_source(pkg)
%STAGE_SOURCE Isolate the smoother's temporary output from the source package.
for k = 1:numel(pkg.source_paths)
    if ~strcmp(surfsmooth3d.edgepreserve.file_sha256(pkg.source_paths{k}),pkg.source_hashes{k})
        error('edgepreserve:sourceChanged', ...
            'Source changed since loading. Browse the package again: %s',pkg.source_paths{k});
    end
end
directory = tempname;
mkdir(directory);
cleanup = onCleanup(@() rmdir(directory,'s'));
scaffold = fullfile(directory,'scaffold.gidmsh');
copyfile(pkg.scaffold_file,scaffold);
source = fullfile(directory,'cad_points.txt');
surfsmooth3d.edgepreserve.write_point_source(pkg.surface,source);
work = struct('directory',directory,'scaffold_file',scaffold,'cad_file',source);
end

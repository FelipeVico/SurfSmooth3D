function export_surface(filename, S, report)
%EXPORT_SURFACE Gate, round-trip check, then publish a .go3 atomically.
if ~report.is_valid
    error('edgepreserve:invalidSurface', ...
        'Refusing invalid geometry: %d bad nodes in %d patches.', ...
        report.stats.bad_node_count,report.stats.bad_patch_count);
end
directory = fileparts(filename);
if isfile(filename)
    error('edgepreserve:existingOutput', 'Refusing to overwrite %s.',filename);
end
temporary = [tempname(directory) '.go3'];
cleanup = onCleanup(@() remove_temporary(temporary));
surfsmooth3d.edgepreserve.write_go3(temporary,S);
reloaded = surfsmooth3d.surfer.load_from_file(temporary);
if reloaded.npts ~= S.npts || reloaded.npatches ~= S.npatches || ...
        ~isequal(reloaded.norders,S.norders) || ~isequal(reloaded.iptype,S.iptype)
    error('edgepreserve:exportLayout', 'Saved surface layout does not match.');
end
[original,~,~,~,~,weights] = extract_arrays(S);
[loaded,~,~,~,~,loadedWeights] = extract_arrays(reloaded);
scale = max(1,max(abs(original(:))));
if any(~isfinite(loaded(:))) || max(abs(original(:)-loaded(:))) > 1e-12*scale || ...
        any(~isfinite(loadedWeights)) || ...
        max(abs(weights(:)-loadedWeights(:))) > 1e-12*max(1,max(abs(weights)))
    error('edgepreserve:exportData', 'Saved surface failed the data round-trip check.');
end
[ok,message] = movefile(temporary,filename);
if ~ok, error('edgepreserve:publish', '%s',message); end
end

function remove_temporary(filename)
if isfile(filename), delete(filename); end
end

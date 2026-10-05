function write_go3(outFile, S)
    [srcvals, ~, ~, ~, ~, ~] = extract_arrays(S);
    norder = infer_common_order(S);

    fid = fopen(outFile, 'w');
    if fid < 0
        error('edgepreserve:go3OpenFailed', ...
            'Could not open output .go3 file for writing: %s', outFile);
    end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '%d\n', norder);
    fprintf(fid, '%d\n', S.npatches);
    fprintf(fid, '%.16e\n', srcvals.');
    clear cleanup;
end


function norder = infer_common_order(S)
    norders = S.norders(:);
    if isempty(norders)
        error('edgepreserve:emptyOrders', ...
            'CAD skeleton has no patch orders.');
    end
    if any(norders ~= norders(1))
        error('edgepreserve:mixedOrders', ...
            'This driver expects one common order across all patches.');
    end
    norder = norders(1);
end

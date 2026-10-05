function write_point_source(S, outFile)
    [srcvals, ~, ~, ~, ~, wts] = extract_arrays(S);
    data = [srcvals(1:3, :).', srcvals(10:12, :).', wts(:)];

    fid = fopen(outFile, 'w');
    if fid < 0
        error('edgepreserve:openFailed', ...
            'Could not open temporary CAD skeleton file for writing: %s', ...
            outFile);
    end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '# x y z nx ny nz w generated from surface.go3\n');
    fprintf(fid, '%.16e %.16e %.16e %.16e %.16e %.16e %.16e\n', data.');
    clear cleanup;
end

function [nodes, tris] = read_gidmsh(fname)
    lines = string(readlines(fname));
    iCoord = find(contains(lines, "Coordinates"), 1, 'first');
    iEndCoord = find(contains(lines, "End Coordinates"), 1, 'first');
    iElem = find(contains(lines, "Elements"), 1, 'first');
    iEndElem = find(contains(lines, "End Elements"), 1, 'first');
    if isempty(iCoord) || isempty(iEndCoord) || isempty(iElem) || ...
            isempty(iEndElem)
        error('edgepreserve:badGidmsh', ...
            'Could not parse GiD mesh sections in %s.', fname);
    end

    nodeLines = lines((iCoord + 1):(iEndCoord - 1));
    elemLines = lines((iElem + 1):(iEndElem - 1));
    nodeData = sscanf(join(nodeLines, newline), '%f');
    elemData = sscanf(join(elemLines, newline), '%f');
    if mod(numel(nodeData), 4) ~= 0 || mod(numel(elemData), 4) ~= 0
        error('edgepreserve:badGidmshData', ...
            'Unexpected GiD data length in %s.', fname);
    end

    nodeData = reshape(nodeData, 4, []);
    elemData = reshape(elemData, 4, []);
    ids = nodeData(1, :);
    if any(ids ~= 1:numel(ids))
        error('edgepreserve:noncontiguousNodes', ...
            'This viewer expects contiguous 1-based GiD node ids.');
    end
    nodes = nodeData(2:4, :);
    tris = elemData(2:4, :);
end

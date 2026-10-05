function [S_blend, beta, dsoft] = blend( ...
        S_cad, S_smooth, edges, params)
    if any(~ismember(params.selectedEdgeIds, [edges.edge_id]))
        error('edgepreserve:unknownEdge', 'Unknown selected CAD edge ID.');
    end
    if isempty(params.selectedEdgeIds)
        dsoft = inf(S_cad.npts, 1);
        beta = zeros(S_cad.npts, 1);
    else
        dsoft = compute_edge_distance(S_cad.r, edges, params);
        beta = beta_from_distance(dsoft, params);
    end

    [srcCad, ~, norders, ~, iptype, ~] = extract_arrays(S_cad);
    [srcSmooth, ~, ~, ~, ~, ~] = extract_arrays(S_smooth);

    betaRow = reshape(beta, 1, []);
    beta3 = repmat(betaRow, 3, 1);
    rBlend = beta3 .* srcSmooth(1:3, :) + ...
        (1 - beta3) .* srcCad(1:3, :);
    srcBlend = surfsmooth3d.edgepreserve.srcvals_from_positions(S_cad, rBlend);

    S_blend = surfsmooth3d.surfer(S_cad.npatches, norders, srcBlend, iptype);
end


function beta = beta_from_distance(dsoft, params)
    sigma = local_sigma(params);
    radius = max(params.c_rad * sigma, eps);
    betaTau = max(params.c_tau * sigma, eps);
    beta = 0.5 * erfc((dsoft(:) - radius) ./ (sqrt(2) * betaTau));
    beta = min(max(beta, 0), 1);
end


function dsoft = compute_edge_distance(xyz, edges, params)
    [edgeXyz, edgeW] = selected_edge_quadrature(edges, ...
        params.selectedEdgeIds);
    if isempty(edgeXyz)
        error('edgepreserve:noSelectedEdges', ...
            'No valid selected CAD edges were found.');
    end

    sigma = local_sigma(params);
    epsSoftmin = params.c_eps * sigma;
    if any(~isfinite(epsSoftmin)) || any(epsSoftmin <= 0)
        error('edgepreserve:badCEps', ...
            'c_eps*sigma must be positive and finite at every node.');
    end

    npts = size(xyz, 2);
    A = zeros(npts, 1);
    chunkSize = 2000;

    for i1 = 1:chunkSize:npts
        i2 = min(npts, i1 + chunkSize - 1);
        xloc = xyz(:, i1:i2).';
        epsLoc = epsSoftmin(i1:i2);
        dx = xloc(:, 1) - edgeXyz(1, :);
        dy = xloc(:, 2) - edgeXyz(2, :);
        dz = xloc(:, 3) - edgeXyz(3, :);
        dist2 = dx.^2 + dy.^2 + dz.^2;
        denom = 2 * epsLoc.^2;
        normalizer = 1 ./ (sqrt(2*pi) * epsLoc);
        A(i1:i2) = normalizer .* ...
            (exp(-dist2 ./ denom) * edgeW(:));
    end

    floorVal = max(params.distanceFloor, realmin);
    dsoft = epsSoftmin .* sqrt(max(0, -2 * log(max(A, floorVal))));
end


function sigma = local_sigma(params)
    sigma = params.sigma(:);
    sigmaFloor = max(params.sigmaFloorRel * params.geometryDiameter, realmin);
    sigma = max(sigma, sigmaFloor);
    if any(~isfinite(sigma)) || any(sigma <= 0)
        error('edgepreserve:badSigma', ...
            'The smoother sigma evaluator returned invalid sigma values.');
    end
end


function [edgeXyz, edgeW] = selected_edge_quadrature(edges, edgeIds)
    edgeXyz = zeros(3, 0);
    edgeW = zeros(1, 0);
    edgeIds = unique(edgeIds(:).');

    for id = edgeIds
        idx = find([edges.edge_id] == id, 1);
        if isempty(idx)
            continue
        end
        panels = edges(idx).panels;
        for j = 1:numel(panels)
            edgeXyz = [edgeXyz, panels(j).xyz]; %#ok<AGROW>
            edgeW = [edgeW, panels(j).w]; %#ok<AGROW>
        end
    end
end

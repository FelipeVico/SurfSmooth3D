function h = plot_surface(S, values, plotArgs, refine)
    if nargin < 2 || isempty(values)
        values = S.mean_curv;
    end
    if nargin < 3 || isempty(plotArgs)
        plotArgs = {};
    end
    if nargin < 4 || isempty(refine)
        refine = 20;
    end

    if any(S.iptype(:) ~= 1)
        h = plot(S, values, plotArgs{:});
        return
    end

    values = expand_surfer_plot_values(S, values);
    fcoefs = S.vals2coefs(values.');
    ncoefs = S.vals2coefs(S.n);
    [uvPlot, triFaces] = local_reference_triangle_mesh(refine);

    npatches = S.npatches;
    nplot = size(uvPlot, 2);
    nfaces = size(triFaces, 1);
    vertices = zeros(npatches*nplot, 3);
    vertexNormals = zeros(npatches*nplot, 3);
    cdata = zeros(npatches*nplot, 1);
    faces = zeros(npatches*nfaces, 3);

    polCache = struct();
    for ipatch = 1:npatches
        order = double(S.norders(ipatch));
        cacheKey = sprintf('o%d', order);
        if ~isfield(polCache, cacheKey)
            polCache.(cacheKey) = surfsmooth3d.internal.koorn.pols(order, uvPlot);
        end
        pols = polCache.(cacheKey);

        i1 = double(S.ixyzs(ipatch));
        i2 = double(S.ixyzs(ipatch + 1)) - 1;
        nodeRange = i1:i2;

        rPlot = S.srccoefs{ipatch}(1:3, :) * pols;
        nPlot = ncoefs(:, nodeRange) * pols;
        nNorm = sqrt(sum(nPlot.^2, 1));
        good = isfinite(nNorm) & nNorm > 0;
        nPlot(:, good) = nPlot(:, good) ./ nNorm(good);
        nPlot(:, ~good) = 0;
        fPlot = fcoefs(nodeRange) * pols;

        vertexRange = (ipatch - 1)*nplot + (1:nplot);
        faceRange = (ipatch - 1)*nfaces + (1:nfaces);
        vertices(vertexRange, :) = rPlot.';
        vertexNormals(vertexRange, :) = nPlot.';
        cdata(vertexRange) = fPlot(:);
        faces(faceRange, :) = triFaces + vertexRange(1) - 1;
    end

    h = patch('Faces', faces, 'Vertices', vertices, ...
        'FaceVertexCData', cdata, ...
        'FaceColor', 'interp', ...
        'EdgeColor', 'none', ...
        'LineStyle', 'none', ...
        'FaceLighting', 'gouraud', ...
        'SpecularStrength', 0.05, ...
        'DiffuseStrength', 0.65, ...
        'AmbientStrength', 0.45, ...
        plotArgs{:});
    set(h, 'VertexNormals', vertexNormals);
    axis equal;
    view(3);
    grid on;
end


function values = expand_surfer_plot_values(S, values)
    values = values(:);
    if numel(values) == S.npts
        return
    end
    if numel(values) == S.npatches
        patchValues = values;
        values = zeros(S.npts, 1);
        for ipatch = 1:S.npatches
            i1 = double(S.ixyzs(ipatch));
            i2 = double(S.ixyzs(ipatch + 1)) - 1;
            values(i1:i2) = patchValues(ipatch);
        end
        return
    end
    error('edgepreserve:badPlotValues', ...
        'Plot values must have length S.npts or S.npatches.');
end


function [uv, faces] = local_reference_triangle_mesh(refine)
    refine = max(1, round(refine));
    index = nan(refine + 1, refine + 1);
    uv = zeros(2, (refine + 1)*(refine + 2)/2);
    k = 0;
    for i = 0:refine
        for j = 0:(refine - i)
            k = k + 1;
            index(i + 1, j + 1) = k;
            uv(:, k) = [i; j] / refine;
        end
    end

    faces = zeros(refine*refine, 3);
    nf = 0;
    for i = 0:(refine - 1)
        for j = 0:(refine - 1 - i)
            v00 = index(i + 1, j + 1);
            v10 = index(i + 2, j + 1);
            v01 = index(i + 1, j + 2);
            nf = nf + 1;
            faces(nf, :) = [v00 v10 v01];
            if j <= refine - 2 - i
                v11 = index(i + 2, j + 2);
                nf = nf + 1;
                faces(nf, :) = [v10 v11 v01];
            end
        end
    end
    faces = faces(1:nf, :);
end

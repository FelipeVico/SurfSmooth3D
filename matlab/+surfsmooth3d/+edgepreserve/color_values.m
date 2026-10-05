function [plotValues, colorLabel, statusNote] = color_values( ...
        S_display, S_cad, beta, colorMode)
    statusNote = '';
    switch colorMode
        case 'beta'
            if numel(beta) == S_display.npts
                plotValues = beta(:);
                colorLabel = 'beta (0 CAD, 1 smooth)';
            else
                plotValues = zeros(S_display.npts, 1);
                colorLabel = 'constant';
                statusNote = 'beta unavailable in fully smooth mode';
            end
        case 'signed mean curvature'
            plotValues = S_display.mean_curv(:);
            colorLabel = 'signed mean curvature';
        case 'log10 abs mean curvature'
            plotValues = log10(abs(S_display.mean_curv(:)) + 1e-16);
            colorLabel = 'log10(|mean curvature| + 1e-16)';
        case 'patch max abs curvature'
            plotValues = patch_max_abs_curvature_values(S_display);
            colorLabel = 'patch max |mean curvature|';
        case 'jacobian ratio to CAD'
            [srcDisplay, ~, ~, ~, ~, ~] = extract_arrays(S_display);
            [srcCad, ~, ~, ~, ~, ~] = extract_arrays(S_cad);
            jacDisplay = surface_jacobian_from_srcvals(srcDisplay);
            jacCad = surface_jacobian_from_srcvals(srcCad);
            plotValues = jacDisplay ./ max(jacCad, realmin);
            colorLabel = '|du x dv| / CAD';
        case 'displacement to CAD'
            plotValues = sqrt(sum((S_display.r - S_cad.r).^2, 1)).';
            colorLabel = 'displacement to CAD';
        case 'constant'
            plotValues = zeros(S_display.npts, 1);
            colorLabel = 'constant';
        otherwise
            error('edgepreserve:badColorMode', ...
                'Unknown color mode "%s".', colorMode);
    end
end


function values = patch_max_abs_curvature_values(S)
    values = accumarray(S.patch_id(:), abs(S.mean_curv(:)), ...
        [S.npatches 1], @max, nan);
end


function jac = surface_jacobian_from_srcvals(srcvals)
    cr = cross(srcvals(4:6, :).', srcvals(7:9, :).').';
    jac = sqrt(sum(cr.^2, 1)).';
end

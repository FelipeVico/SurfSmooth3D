function limitText = color_limits(ax, plotValues, colorMode, statusNote)
    if strcmp(colorMode, 'signed mean curvature')
        finiteVals = plotValues(isfinite(plotValues));
        if isempty(finiteVals)
            climVals = [-1 1];
            statsText = 'signed H unavailable';
        else
            maxAbs = max(abs(finiteVals));
            if maxAbs <= 0
                maxAbs = 1;
            end
            climVals = [-maxAbs maxAbs];
            statsText = sprintf('signed H min/mean/max %.3g %.3g %.3g', ...
                min(finiteVals), mean(finiteVals), max(finiteVals));
        end
        clim(ax, climVals);
        limitText = sprintf('\n%s | clim [%.3g %.3g]', ...
            statsText, climVals(1), climVals(2));
        return
    end

    if strcmp(colorMode, 'beta') && isempty(statusNote)
        clim(ax, [0 1]);
        limitText = sprintf('\nclim [0 1]');
        return
    end

    if strcmp(colorMode, 'constant') || ...
            strcmp(statusNote, 'beta unavailable in fully smooth mode')
        clim(ax, [0 1]);
        limitText = sprintf('\nclim [0 1]');
        return
    end

    finiteVals = plotValues(isfinite(plotValues));
    if isempty(finiteVals)
        clim(ax, [0 1]);
        limitText = sprintf('\nclim [0 1]');
        return
    end

    lo = min(finiteVals);
    hi = max(finiteVals);
    if lo == hi
        pad = max(1, abs(lo))*1e-6;
        lo = lo - pad;
        hi = hi + pad;
    end
    clim(ax, [lo hi]);
    limitText = sprintf('\nclim [%.3g %.3g]', lo, hi);
end

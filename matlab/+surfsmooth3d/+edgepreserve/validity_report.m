function report = validity_report(S, label, S_ref)
    if nargin < 3
        S_ref = [];
    end

    thresholds = struct();
    thresholds.normal_cross_dot_min = 0.99;
    thresholds.jac_ratio_fail_min = 1e-3;
    thresholds.jac_ratio_warn_min = 1e-2;

    [srcvals, ~, ~, ~, ~, wts] = extract_arrays(S);
    cr = cross(srcvals(4:6, :).', srcvals(7:9, :).').';
    jac = sqrt(sum(cr.^2, 1)).';
    crUnit = cr ./ max(jac.', realmin);
    normalCrossDot = sum(srcvals(10:12, :) .* crUnit, 1).';

    hasRef = ~isempty(S_ref) && S_ref.npts == S.npts && ...
        S_ref.npatches == S.npatches;
    if hasRef
        [srcRef, ~, ~, ~, ~, ~] = extract_arrays(S_ref);
        crRef = cross(srcRef(4:6, :).', srcRef(7:9, :).').';
        jacRef = sqrt(sum(crRef.^2, 1)).';
        jacRatioToRef = jac ./ max(jacRef, realmin);
        displacementToRef = sqrt(sum((S.r - S_ref.r).^2, 1)).';
        curvatureThreshold = 100 * max(1, max(abs(S_ref.mean_curv(:))));
    else
        jacRatioToRef = nan(S.npts, 1);
        displacementToRef = nan(S.npts, 1);
        curvatureThreshold = inf;
    end
    thresholds.curvature_abs_max = curvatureThreshold;

    finiteNode = all(isfinite(srcvals), 1).' & isfinite(wts(:)) & ...
        isfinite(S.mean_curv(:)) & isfinite(jac) & ...
        isfinite(normalCrossDot);
    if hasRef
        finiteNode = finiteNode & isfinite(jacRatioToRef) & ...
            isfinite(displacementToRef);
    end

    badFinite = ~finiteNode;
    badNormal = normalCrossDot < thresholds.normal_cross_dot_min;
    badJacRatio = hasRef & ...
        (jacRatioToRef < thresholds.jac_ratio_fail_min);
    warnJacRatio = hasRef & ...
        (jacRatioToRef < thresholds.jac_ratio_warn_min) & ...
        ~badJacRatio;
    badCurvature = abs(S.mean_curv(:)) > thresholds.curvature_abs_max;
    badMask = badFinite | badNormal | badJacRatio | badCurvature;

    failureCode = zeros(S.npts, 1);
    failureCode(warnJacRatio) = 1;
    failureCode(badCurvature) = 2;
    failureCode(badJacRatio) = 3;
    failureCode(badNormal) = 4;
    failureCode(badFinite) = 5;

    report = struct();
    report.label = char(label);
    report.thresholds = thresholds;
    report.has_reference = hasRef;
    report.normal_cross_dot = normalCrossDot;
    report.jac = jac;
    report.jac_ratio_to_ref = jacRatioToRef;
    report.displacement_to_ref = displacementToRef;
    report.r = S.r;
    report.patch_id = S.patch_id(:);
    report.wts = wts(:);
    report.mean_curv = S.mean_curv(:);
    report.bad_mask = badMask;
    report.warning_mask = warnJacRatio;
    report.failure_code = failureCode;
    report.is_valid = ~any(badMask);
    report.has_warnings = any(warnJacRatio);
    report.csv_file = fullfile(tempdir,'surfsmooth3d_surface_validity_bad_nodes.csv');

    report.stats = struct();
    report.stats.npts = S.npts;
    report.stats.npatches = S.npatches;
    report.stats.normal_cross_dot_min = min(normalCrossDot);
    report.stats.normal_cross_dot_mean = mean(normalCrossDot);
    report.stats.normal_cross_dot_max = max(normalCrossDot);
    report.stats.jac_min = min(jac);
    report.stats.jac_mean = mean(jac);
    report.stats.jac_max = max(jac);
    report.stats.wts_min = min(wts);
    report.stats.wts_mean = mean(wts);
    report.stats.wts_max = max(wts);
    report.stats.mean_curv_min = min(S.mean_curv);
    report.stats.mean_curv_mean = mean(S.mean_curv);
    report.stats.mean_curv_max = max(S.mean_curv);
    report.stats.mean_curv_std = std(S.mean_curv);
    report.stats.mean_curv_abs_max = max(abs(S.mean_curv));
    report.stats.bad_node_count = nnz(badMask);
    report.stats.warning_node_count = nnz(warnJacRatio);
    report.stats.bad_patch_count = count_patches(S, badMask);
    report.stats.warning_patch_count = count_patches(S, warnJacRatio);
    report.stats.bad_finite_count = nnz(badFinite);
    report.stats.bad_normal_count = nnz(badNormal);
    report.stats.bad_jac_ratio_count = nnz(badJacRatio);
    report.stats.warn_jac_ratio_count = nnz(warnJacRatio);
    report.stats.bad_curvature_count = nnz(badCurvature);
    if hasRef
        report.stats.jac_ratio_min = min(jacRatioToRef);
        report.stats.jac_ratio_mean = mean(jacRatioToRef);
        report.stats.jac_ratio_max = max(jacRatioToRef);
        report.stats.displacement_min = min(displacementToRef);
        report.stats.displacement_mean = mean(displacementToRef);
        report.stats.displacement_max = max(displacementToRef);
        report.stats.displacement_std = std(displacementToRef);
    else
        report.stats.jac_ratio_min = nan;
        report.stats.jac_ratio_mean = nan;
        report.stats.jac_ratio_max = nan;
        report.stats.displacement_min = nan;
        report.stats.displacement_mean = nan;
        report.stats.displacement_max = nan;
        report.stats.displacement_std = nan;
    end
end


function n = count_patches(S, mask)
    if any(mask)
        n = numel(unique(S.patch_id(mask)));
    else
        n = 0;
    end
end

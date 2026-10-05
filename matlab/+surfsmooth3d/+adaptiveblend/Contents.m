% Adaptive refinement of partially smoothed CAD surfaces.
%
% launch_gui       - Independent GUI; selected CAD edges are smoothed.
% default_settings - Fixed order, sigma and blend controls, tolerance and limits.
% solve            - Two-stage initialization and selective four-way refinement.
% save_result      - Verified partial/smooth pair export and provenance.
%
% The tolerance applies to the partial surface, not its smooth companion.
% Coefficient tails and independent position/normal checks use the fixed CAD
% polynomials, original-scaffold sigma field, and unchanged edge quadrature.
% Independent checks project along inherited launch lines, not along normals
% of the blended output. Their normal reference includes derivatives of beta.
%
% Geometric degeneracy and nondifferentiable soft-distance clamp joins fail
% checks rather than receiving zero error. A depth/point limit returns the
% last completed mesh with converged=false; valid exports get _tolunmet.
% The existing blend formula itself can limit convergence: no formula change
% or automatic regularization is performed here.
%
% Hanging reference edges are allowed. Exact interpatch continuity is not
% enforced, and the tolerance is not a certified BIE solution-error bound.
% Native sessions exist only during solve and are explicitly closed by
% onCleanup, including after a reported Newton failure.

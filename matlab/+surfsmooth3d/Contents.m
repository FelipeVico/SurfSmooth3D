% SurfSmooth3D: CAD geometry, surface smoothing, and .go3 export.
%
% Setup: addpath('matlab'); setup_surfsmooth3d
%
% Smoother interfaces (the original arguments/options are retained):
%   surfsmooth3d.multiscale_mesher          - One/two-stage uniform smoothing.
%   surfsmooth3d.multiscale_mesher_adaptive - Adaptive two-stage smoothing.
%   surfsmooth3d.multiscale_mesher_sigma_eval - Sigma and its spatial gradient.
%   surfsmooth3d.surfer                    - High-order surface representation.
%
% CAD and workflow packages:
%   surfsmooth3d.stepmesher      - STEP meshing and geometry package I/O.
%   surfsmooth3d.edgepreserve    - Fixed-CAD blending, preview, and batch export.
%   surfsmooth3d.adaptivesmoother - Adaptive fully smoothed surfaces.
%   surfsmooth3d.adaptiveblend   - Adaptive partially smoothed surfaces.
%
% Interactive drivers are in examples/matlab; STEP drivers in examples/step.
% Use surfsmooth3d.root() to locate the installed package directory.
% Internal polynomial helpers are implementation details, not a public API.

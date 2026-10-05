% Fixed-CAD two-stage smoothing preview and batch exports.
%
% Start: examples/matlab/driver_two_stage_edge_smooth_batch.m
% Select ONE package folder containing surface.go3, scaffold.gidmsh,
% edges_gll.txt and metadata.txt. The source stays fixed across all outputs.
% Select CAD edge IDs to SMOOTH; empty IDs give a CAD-only partial surface.
% Tune rlam/adapt_sigma and preview resolution, then Draw / Recompute.
% Edge IDs and c_eps/c_rad/c_tau update only the partial preview live.
% Save Current Pair saves the latest drawn full and partial surfaces, not
% pending expensive settings. Generate And Save Batch snapshots the controls.
% Each run gets parameters.mat, batch_manifest.csv, and valid .go3 outputs.
% Invalid full/partial candidates are rejected independently and diagnosed.
%
% Numerical source conventions:
%   STEP: lattice subdivision, consecutive children for each coarse parent.
%   Smoother: recursive children (u/2,v/2), (.5-u/2,.5-v/2),
%             (.5+u/2,v/2), (u/2,.5+v/2).
% CAD values are matched through reference coordinates, never nearest points.
% nquad is the fixed source order. Sigma uses the original coarse scaffold.
% Changing output resolution therefore does not redefine the level-set source.
% The fixed source's quadrature/approximation error still limits convergence.
% default_input_root finds outputs/step and the included CAD test fixtures.
% Other package folders can be selected.
%
% Tests (after matlab/setup_surfsmooth3d.m and addpath('tests/matlab')):
%   test_edgepreserve_fixed_cad
%   smoke_edgepreserve_real_package
%   smoke_edgepreserve_gui_batch
% The latter performs real Fortran solves and saves test outputs under tempdir.

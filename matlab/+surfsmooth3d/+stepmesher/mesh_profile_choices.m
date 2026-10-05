function choices = mesh_profile_choices()
%MESH_PROFILE_CHOICES Return documented Gmsh scaffold meshing presets.
%
% The "name" column is the user-facing launch choice used by the example
% drivers. The "mesh_profile" column is the lower-level C++/Gmsh profile.

names = [
    "rigid"
    "quality"
    "adaptive"
    "robust"
    "curvature_coarse"
    "curvature_balanced"
    "curvature_balanced_local"
    "curvature_balanced_local_cadedge"
    "curvature_balanced_wide"
    "curvature_fine"
    "curvature_guarded"
    ];

descriptions = [
    "uniform size, original rigid/default profile"
    "uniform size, frontal-Delaunay Gmsh surface algorithm"
    "uniform size, Delaunay Gmsh surface algorithm"
    "uniform size, automatic/robust Gmsh surface algorithm"
    "curvature sizing with relaxed min/max size bounds"
    "curvature sizing, recommended first adaptive test"
    "curvature sizing like balanced, without boundary size spreading"
    "curvature sizing with hmin from one quarter of the smallest meaningful CAD edge"
    "curvature sizing like balanced, but with larger flat-region cap"
    "curvature sizing with smaller hmin/hmax for small features"
    "curvature sizing plus quality enforcement"
    ];

meshProfiles = [
    "rigid"
    "quality"
    "adaptive"
    "robust"
    "curvature"
    "curvature"
    "curvature"
    "curvature"
    "curvature"
    "curvature"
    "curvature"
    ];

curvaturePoints = [50; 50; 50; 50; 35; 50; 50; 50; 50; 80; 50];
hminModes = [
    "fraction"
    "fraction"
    "fraction"
    "fraction"
    "fraction"
    "fraction"
    "fraction"
    "cad_edge"
    "fraction"
    "fraction"
    "fraction"
    ];
hminFractions = [0.01; 0.01; 0.01; 0.01; 0.02; 0.01; 0.01; 0.01; 0.01; 0.005; 0.01];
hminEdgeScales = [0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25];
cadEdgeTolFractions = [1e-6; 1e-6; 1e-6; 1e-6; 1e-6; 1e-6; 1e-6; 1e-6; 1e-6; 1e-6; 1e-6];
hmaxFractions = [0.10; 0.10; 0.10; 0.10; 0.15; 0.10; 0.10; 0.10; 0.20; 0.05; 0.10];
meshSizeExtendFromBoundary = [
    true
    true
    true
    true
    true
    true
    false
    false
    true
    true
    true
    ];
gmshAlgorithms = [
    "frontal_delaunay"
    "frontal_delaunay"
    "delaunay"
    "automatic"
    "frontal_delaunay"
    "frontal_delaunay"
    "frontal_delaunay"
    "frontal_delaunay"
    "frontal_delaunay"
    "frontal_delaunay"
    "frontal_delaunay"
    ];
optimize = [true; true; true; true; true; true; true; true; true; true; true];
minAngleDeg = [15; 15; 15; 15; 12; 15; 15; 15; 15; 18; 15];
maxEdgeRatio = [6; 6; 6; 6; 8; 6; 6; 6; 6; 5; 6];
enforceQuality = [false; false; false; false; false; false; false; false; false; false; true];

choices = table(names, descriptions, meshProfiles, curvaturePoints, ...
    hminModes, hminFractions, hminEdgeScales, cadEdgeTolFractions, ...
    hmaxFractions, meshSizeExtendFromBoundary, ...
    gmshAlgorithms, optimize, minAngleDeg, maxEdgeRatio, enforceQuality, ...
    'VariableNames', { ...
    'name', 'description', 'mesh_profile', 'curvature_points', ...
    'hmin_mode', 'hmin_fraction', 'hmin_edge_scale', ...
    'cad_edge_tol_fraction', 'hmax_fraction', ...
    'mesh_size_extend_from_boundary', ...
    'gmsh_algorithm', 'optimize', 'min_angle_deg', 'max_edge_ratio', ...
    'enforce_quality'});
end

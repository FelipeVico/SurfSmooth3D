#pragma once

#include <string>

namespace stepmesher {

struct MeshOptions {
  int order = 8;
  double mesh_fraction = 0.10;
  int refinement_level = 0;
  int edge_gll_order = 16;
  std::string mesh_profile = "rigid";
  int curvature_points = 50;
  std::string hmin_mode = "fraction";
  double hmin_fraction = 0.01;
  double hmin_edge_scale = 0.25;
  double cad_edge_tol_fraction = 1.0e-6;
  double hmax_fraction = 0.10;
  bool mesh_size_extend_from_boundary = true;
  std::string gmsh_algorithm = "frontal_delaunay";
  bool optimize = true;
  double min_angle_deg = 15.0;
  double max_edge_ratio = 6.0;
  bool enforce_quality = false;
  double occt_precision = 1.0e-6;
  double occt_maxprecision = 1.0e-6;
  std::string occ_info_exe;
  bool sameparameter = true;
  std::string surfacecurve_mode = "3d_preferred";
  bool verbose = false;
};

} // namespace stepmesher

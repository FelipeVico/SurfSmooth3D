#pragma once

#include <string>
#include <vector>

namespace stepmesher {

struct MeshResult {
  int native_identity_faces = 0;
  int native_identity_occurrences = 0;
  std::string cad_identity = "native_occ_topology_v1";
  std::string occ_evaluation;
  std::string runtime_version = "3.0.0";
  std::string gmsh_version = "4.12.2";
  std::string occt_version = "7.9.3";
  int order = 0;
  int iptype = 1;
  int nodes_per_patch = 0;
  int npatches = 0;
  int nscaffold_nodes = 0;
  int nscaffold_triangles = 0;
  int projection_failures = 0;
  int step1b_edge_newton_success = 0;
  int step1b_edge_curve_projection_success = 0;
  int edge_newton_success = 0;
  int edge_curve_projection_success = 0;
  int edge_fallback_points = 0;
  int continuity_break_count = 0;
  int continuity_helper_curve_count = 0;
  int fragmented_face_count = 0;
  int feature_edge_fallback_points = 0;
  int edge_order = 0;
  int edge_nodes_per_panel = 0;
  int ncad_edges = 0;
  int npanels_total = 0;
  int bad_min_angle_count = 0;
  int bad_edge_ratio_count = 0;
  int worst_quality_triangle = 0;

  double bbox_diag = 0.0;
  double mesh_size = 0.0;
  double hmin_effective = 0.0;
  double hmax_effective = 0.0;
  double cad_area = 0.0;
  double min_raw_cad_edge_length = 0.0;
  double min_meaningful_cad_edge_length = 0.0;
  double max_projection_distance = 0.0;
  double mean_projection_distance = 0.0;
  double min_scaffold_angle_deg = 0.0;
  double max_scaffold_edge_ratio = 0.0;
  int cad_edge_count = 0;
  int meaningful_cad_edge_count = 0;

  std::string step_file;
  std::string mesh_profile;
  std::string sample_rule;
  std::string backend = "gmsh_occ";

  // MATLAB-facing arrays are stored in column-major order.
  // xyz has shape 3 x (nodes_per_patch * npatches).
  std::vector<double> xyz;
  std::vector<double> proj_dist;

  // One entry per high-order patch.
  std::vector<int> tri_face_tags;
  std::vector<int> parent_coarse_tri;

  // Scaffold/debug data, also MATLAB-facing and 1-based where relevant.
  std::vector<double> scaffold_nodes;    // 3 x nscaffold_nodes
  std::vector<int> scaffold_node_tags;   // original Gmsh node tags
  std::vector<int> scaffold_triangles;   // local node ids, 3 x nscaffold_triangles
  std::vector<int> scaffold_triangle_node_tags; // original Gmsh tags
  std::vector<int> scaffold_face_tags;   // 1 x nscaffold_triangles
  std::vector<double> scaffold_min_angle_deg; // 1 x nscaffold_triangles
  std::vector<double> scaffold_edge_ratio;    // 1 x nscaffold_triangles
  std::vector<double> scaffold_triangle_area; // 1 x nscaffold_triangles

  // High-order CAD edges with piecewise GLL panels, also MATLAB-facing.
  std::vector<double> edge_gll_nodes;     // 1 x edge_nodes_per_panel
  std::vector<int> edge_ids;              // 1 x ncad_edges
  std::vector<int> edge_curve_tags;       // 1 x ncad_edges
  std::vector<int> edge_adjacent_faces;   // 2 x ncad_edges, -1 when unknown
  std::vector<int> edge_panel_start;      // 1 x ncad_edges, 1-based panel index
  std::vector<int> edge_panel_count;      // 1 x ncad_edges
  std::vector<int> panel_node_tags;       // 2 x npanels_total
  std::vector<double> panel_xyz;          // 3 x (edge_nodes_per_panel*npanels_total)
  std::vector<double> panel_tangents;     // 3 x (edge_nodes_per_panel*npanels_total)
  std::vector<double> panel_weights;      // 1 x (edge_nodes_per_panel*npanels_total)
};

} // namespace stepmesher

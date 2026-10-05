#include "step_mesher/mesher.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void usage() {
  std::cerr
      << "Usage:\n"
      << "  step_mesher_cli mesh --step file.step --order 8 "
         "--mesh-fraction 0.10 --refinement 0 --edge-order 16 "
         "--profile rigid --occ-info-exe helper --out outdir\n"
      << "  Curvature options include --curvature-points, --hmin-fraction, "
         "--hmin-mode, --hmin-edge-scale, --cad-edge-tol-fraction, "
         "--hmax-fraction, and --no-mesh-size-extend-from-boundary.\n";
}

std::string next_arg(int &i, int argc, char **argv, const std::string &name) {
  if (i + 1 >= argc) {
    throw std::runtime_error("Missing value for " + name);
  }
  ++i;
  return argv[i];
}

void write_meta(const std::filesystem::path &path,
                const stepmesher::MeshResult &mesh,
                const stepmesher::MeshOptions &opts) {
  std::ofstream f(path);
  f << std::setprecision(17);
  f << "{\n";
  f << "  \"cad_identity\": \"" << mesh.cad_identity << "\",\n";
  f << "  \"occ_evaluation\": \"" << mesh.occ_evaluation << "\",\n";
  f << "  \"runtime_version\": \"" << mesh.runtime_version << "\",\n";
  f << "  \"gmsh_version\": \"" << mesh.gmsh_version << "\",\n";
  f << "  \"occt_version\": \"" << mesh.occt_version << "\",\n";
  f << "  \"native_identity_faces\": " << mesh.native_identity_faces << ",\n";
  f << "  \"native_identity_occurrences\": " << mesh.native_identity_occurrences << ",\n";

  f << "  \"backend\": \"" << mesh.backend << "\",\n";
  f << "  \"step_file\": \"" << mesh.step_file << "\",\n";
  f << "  \"order\": " << mesh.order << ",\n";
  f << "  \"iptype\": " << mesh.iptype << ",\n";
  f << "  \"mesh_fraction\": " << opts.mesh_fraction << ",\n";
  f << "  \"refinement_level\": " << opts.refinement_level << ",\n";
  f << "  \"edge_gll_order\": " << opts.edge_gll_order << ",\n";
  f << "  \"mesh_profile\": \"" << opts.mesh_profile << "\",\n";
  f << "  \"curvature_points\": " << opts.curvature_points << ",\n";
  f << "  \"hmin_mode\": \"" << opts.hmin_mode << "\",\n";
  f << "  \"hmin_fraction\": " << opts.hmin_fraction << ",\n";
  f << "  \"hmin_edge_scale\": " << opts.hmin_edge_scale << ",\n";
  f << "  \"cad_edge_tol_fraction\": " << opts.cad_edge_tol_fraction << ",\n";
  f << "  \"hmax_fraction\": " << opts.hmax_fraction << ",\n";
  f << "  \"hmin_effective\": " << mesh.hmin_effective << ",\n";
  f << "  \"hmax_effective\": " << mesh.hmax_effective << ",\n";
  f << "  \"mesh_size_extend_from_boundary\": "
    << (opts.mesh_size_extend_from_boundary ? "true" : "false") << ",\n";
  f << "  \"gmsh_algorithm\": \"" << opts.gmsh_algorithm << "\",\n";
  f << "  \"optimize\": " << (opts.optimize ? "true" : "false") << ",\n";
  f << "  \"min_angle_deg\": " << opts.min_angle_deg << ",\n";
  f << "  \"max_edge_ratio\": " << opts.max_edge_ratio << ",\n";
  f << "  \"enforce_quality\": "
    << (opts.enforce_quality ? "true" : "false") << ",\n";
  f << "  \"occ_info_exe\": \"" << opts.occ_info_exe << "\",\n";
  f << "  \"sample_rule\": \"equispaced_cli_fallback\",\n";
  f << "  \"nodes_per_patch\": " << mesh.nodes_per_patch << ",\n";
  f << "  \"npatches\": " << mesh.npatches << ",\n";
  f << "  \"nscaffold_nodes\": " << mesh.nscaffold_nodes << ",\n";
  f << "  \"nscaffold_triangles\": " << mesh.nscaffold_triangles << ",\n";
  f << "  \"bbox_diag\": " << mesh.bbox_diag << ",\n";
  f << "  \"mesh_size\": " << mesh.mesh_size << ",\n";
  f << "  \"min_raw_cad_edge_length\": " << mesh.min_raw_cad_edge_length
    << ",\n";
  f << "  \"min_meaningful_cad_edge_length\": "
    << mesh.min_meaningful_cad_edge_length << ",\n";
  f << "  \"cad_edge_count\": " << mesh.cad_edge_count << ",\n";
  f << "  \"meaningful_cad_edge_count\": " << mesh.meaningful_cad_edge_count
    << ",\n";
  f << "  \"cad_area\": " << mesh.cad_area << ",\n";
  f << "  \"max_projection_distance\": " << mesh.max_projection_distance
    << ",\n";
  f << "  \"mean_projection_distance\": " << mesh.mean_projection_distance
    << ",\n";
  f << "  \"projection_failures\": " << mesh.projection_failures << ",\n";
  f << "  \"min_scaffold_angle_deg\": " << mesh.min_scaffold_angle_deg
    << ",\n";
  f << "  \"max_scaffold_edge_ratio\": " << mesh.max_scaffold_edge_ratio
    << ",\n";
  f << "  \"bad_min_angle_count\": " << mesh.bad_min_angle_count << ",\n";
  f << "  \"bad_edge_ratio_count\": " << mesh.bad_edge_ratio_count << ",\n";
  f << "  \"worst_quality_triangle\": " << mesh.worst_quality_triangle
    << ",\n";
  f << "  \"step1b_edge_newton_success\": "
    << mesh.step1b_edge_newton_success << ",\n";
  f << "  \"step1b_edge_curve_projection_success\": "
    << mesh.step1b_edge_curve_projection_success << ",\n";
  f << "  \"edge_newton_success\": " << mesh.edge_newton_success << ",\n";
  f << "  \"edge_curve_projection_success\": "
    << mesh.edge_curve_projection_success << ",\n";
  f << "  \"edge_fallback_points\": " << mesh.edge_fallback_points << ",\n";
  f << "  \"continuity_break_count\": " << mesh.continuity_break_count
    << ",\n";
  f << "  \"continuity_helper_curve_count\": "
    << mesh.continuity_helper_curve_count << ",\n";
  f << "  \"fragmented_face_count\": " << mesh.fragmented_face_count << ",\n";
  f << "  \"feature_edge_fallback_points\": "
    << mesh.feature_edge_fallback_points << "\n";
  f << "}\n";
}

void write_xyz_csv(const std::filesystem::path &path,
                   const stepmesher::MeshResult &mesh) {
  std::ofstream f(path);
  f << std::setprecision(17);
  f << "patch,node,x,y,z,proj_dist,face_tag,parent_coarse_tri\n";
  for (int patch = 0; patch < mesh.npatches; ++patch) {
    for (int node = 0; node < mesh.nodes_per_patch; ++node) {
      const std::size_t col =
          static_cast<std::size_t>(patch * mesh.nodes_per_patch + node);
      f << patch + 1 << "," << node + 1 << "," << mesh.xyz[3 * col] << ","
        << mesh.xyz[3 * col + 1] << "," << mesh.xyz[3 * col + 2] << ","
        << mesh.proj_dist[col] << "," << mesh.tri_face_tags[patch] << ","
        << mesh.parent_coarse_tri[patch] << "\n";
    }
  }
}

void write_scaffold_csv(const std::filesystem::path &nodes_path,
                        const std::filesystem::path &tri_path,
                        const stepmesher::MeshResult &mesh) {
  std::ofstream nf(nodes_path);
  nf << std::setprecision(17);
  nf << "node,gmsh_tag,x,y,z\n";
  for (int i = 0; i < mesh.nscaffold_nodes; ++i) {
    nf << i + 1 << "," << mesh.scaffold_node_tags[i] << ","
       << mesh.scaffold_nodes[3 * i] << ","
       << mesh.scaffold_nodes[3 * i + 1] << ","
       << mesh.scaffold_nodes[3 * i + 2] << "\n";
  }

  std::ofstream tf(tri_path);
  tf << "triangle,node1,node2,node3,node1_tag,node2_tag,node3_tag,"
        "face_tag,min_angle_deg,edge_ratio,area\n";
  for (int i = 0; i < mesh.nscaffold_triangles; ++i) {
    tf << i + 1 << "," << mesh.scaffold_triangles[3 * i] << ","
       << mesh.scaffold_triangles[3 * i + 1] << ","
       << mesh.scaffold_triangles[3 * i + 2] << ","
       << mesh.scaffold_triangle_node_tags[3 * i] << ","
       << mesh.scaffold_triangle_node_tags[3 * i + 1] << ","
       << mesh.scaffold_triangle_node_tags[3 * i + 2] << ","
       << mesh.scaffold_face_tags[i] << ","
       << mesh.scaffold_min_angle_deg[i] << ","
       << mesh.scaffold_edge_ratio[i] << ","
       << mesh.scaffold_triangle_area[i] << "\n";
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2 || std::string(argv[1]) != "mesh") {
    usage();
    return 2;
  }

  stepmesher::MeshOptions opts;
  std::string step_file;
  std::filesystem::path outdir = "step_mesher_out";

  try {
    for (int i = 2; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--step") {
        step_file = next_arg(i, argc, argv, arg);
      } else if (arg == "--order") {
        opts.order = std::stoi(next_arg(i, argc, argv, arg));
      } else if (arg == "--mesh-fraction") {
        opts.mesh_fraction = std::stod(next_arg(i, argc, argv, arg));
      } else if (arg == "--refinement") {
        opts.refinement_level = std::stoi(next_arg(i, argc, argv, arg));
      } else if (arg == "--edge-order") {
        opts.edge_gll_order = std::stoi(next_arg(i, argc, argv, arg));
      } else if (arg == "--profile") {
        opts.mesh_profile = next_arg(i, argc, argv, arg);
      } else if (arg == "--curvature-points") {
        opts.curvature_points = std::stoi(next_arg(i, argc, argv, arg));
      } else if (arg == "--hmin-mode") {
        opts.hmin_mode = next_arg(i, argc, argv, arg);
      } else if (arg == "--hmin-fraction") {
        opts.hmin_fraction = std::stod(next_arg(i, argc, argv, arg));
      } else if (arg == "--hmin-edge-scale") {
        opts.hmin_edge_scale = std::stod(next_arg(i, argc, argv, arg));
      } else if (arg == "--cad-edge-tol-fraction") {
        opts.cad_edge_tol_fraction = std::stod(next_arg(i, argc, argv, arg));
      } else if (arg == "--hmax-fraction") {
        opts.hmax_fraction = std::stod(next_arg(i, argc, argv, arg));
      } else if (arg == "--mesh-size-extend-from-boundary") {
        opts.mesh_size_extend_from_boundary =
            std::stod(next_arg(i, argc, argv, arg)) != 0.0;
      } else if (arg == "--no-mesh-size-extend-from-boundary") {
        opts.mesh_size_extend_from_boundary = false;
      } else if (arg == "--gmsh-algorithm") {
        opts.gmsh_algorithm = next_arg(i, argc, argv, arg);
      } else if (arg == "--no-optimize") {
        opts.optimize = false;
      } else if (arg == "--min-angle") {
        opts.min_angle_deg = std::stod(next_arg(i, argc, argv, arg));
      } else if (arg == "--max-edge-ratio") {
        opts.max_edge_ratio = std::stod(next_arg(i, argc, argv, arg));
      } else if (arg == "--enforce-quality") {
        opts.enforce_quality = true;
      } else if (arg == "--occ-info-exe") {
        opts.occ_info_exe = next_arg(i, argc, argv, arg);
      } else if (arg == "--out") {
        outdir = next_arg(i, argc, argv, arg);
      } else if (arg == "--verbose") {
        opts.verbose = true;
      } else {
        throw std::runtime_error("Unknown argument: " + arg);
      }
    }

    if (step_file.empty()) {
      throw std::runtime_error("--step is required");
    }

    std::filesystem::create_directories(outdir);
    auto samples = stepmesher::equispaced_triangle_nodes(opts.order);
    auto mesh = stepmesher::mesh_step_file(step_file, opts, samples);
    mesh.sample_rule = "equispaced_cli_fallback";

    write_meta(outdir / "meta.json", mesh, opts);
    write_xyz_csv(outdir / "xyz.csv", mesh);
    write_scaffold_csv(outdir / "scaffold_nodes.csv",
                       outdir / "scaffold_triangles.csv", mesh);

    std::cout << "Wrote " << mesh.npatches << " patches to " << outdir
              << "\n";
    std::cout << "Note: CLI sampling uses an equispaced fallback. Use MATLAB "
                 "export_xyzw for RV quadrature xyzw output.\n";
  } catch (const std::exception &e) {
    std::cerr << "step_mesher_cli: " << e.what() << "\n";
    return 1;
  }

  return 0;
}

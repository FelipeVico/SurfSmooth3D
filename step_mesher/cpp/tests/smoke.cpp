#include "step_mesher/mesher.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: step_mesher_smoke file.step\n";
    return 2;
  }

  try {
    stepmesher::MeshOptions opts;
    opts.order = 2;
    opts.mesh_fraction = 0.20;
    opts.refinement_level = 0;
    auto samples = stepmesher::equispaced_triangle_nodes(opts.order);
    auto mesh = stepmesher::mesh_step_file(argv[1], opts, samples);
    if (mesh.npatches <= 0 || mesh.xyz.empty()) {
      throw std::runtime_error("Mesh result is empty");
    }
    std::cout << "patches=" << mesh.npatches
              << " nodes_per_patch=" << mesh.nodes_per_patch
              << " cad_area=" << mesh.cad_area << "\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }

  return 0;
}

#pragma once

#include <vector>

namespace stepmesher {

struct SampleNodes {
  int order = 0;
  int nodes_per_patch = 0;
  std::vector<double> uv; // 2 x nodes_per_patch, column-major
};

int triangle_node_count(int order);
SampleNodes equispaced_triangle_nodes(int order);

} // namespace stepmesher

#include "step_mesher/sample_nodes.hpp"

#include <stdexcept>

namespace stepmesher {

int triangle_node_count(int order) {
  if (order < 0) {
    throw std::invalid_argument("Triangle order must be nonnegative");
  }
  return (order + 1) * (order + 2) / 2;
}

SampleNodes equispaced_triangle_nodes(int order) {
  SampleNodes nodes;
  nodes.order = order;
  nodes.nodes_per_patch = triangle_node_count(order);
  nodes.uv.resize(2 * static_cast<std::size_t>(nodes.nodes_per_patch));

  if (order == 0) {
    nodes.uv[0] = 1.0 / 3.0;
    nodes.uv[1] = 1.0 / 3.0;
    return nodes;
  }

  int k = 0;
  for (int j = 0; j <= order; ++j) {
    for (int i = 0; i <= order - j; ++i) {
      nodes.uv[2 * static_cast<std::size_t>(k)] =
          static_cast<double>(i) / static_cast<double>(order);
      nodes.uv[2 * static_cast<std::size_t>(k) + 1] =
          static_cast<double>(j) / static_cast<double>(order);
      ++k;
    }
  }
  return nodes;
}

} // namespace stepmesher

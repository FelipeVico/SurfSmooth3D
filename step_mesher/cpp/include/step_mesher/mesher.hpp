#pragma once

#include "step_mesher/mesh_result.hpp"
#include "step_mesher/options.hpp"
#include "step_mesher/sample_nodes.hpp"

#include <string>

namespace stepmesher {

MeshResult mesh_step_file(const std::string &step_file,
                          const MeshOptions &options,
                          const SampleNodes &sample_nodes);

} // namespace stepmesher

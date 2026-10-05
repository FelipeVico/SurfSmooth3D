#include "step_mesher/mesher.hpp"
#include "../src/occ_operation.hpp"
#include <iostream>
#include <stdexcept>
#include <filesystem>

void require(bool condition, const char *what) {
  if (!condition) throw std::runtime_error(what);
}

int main(int argc, char **argv) {
  try {
    using stepmesher::detail::call_legacy_occ_operation;
    const auto success = call_legacy_occ_operation([] { return 17; });
    require(success && *success == 17, "operation boundary changed success");
    require(!call_legacy_occ_operation([]() -> int {
      throw std::runtime_error("synthetic helper operation failure");
    }), "operation exception did not preserve legacy fallback");
    require(!call_legacy_occ_operation([]() -> int { throw 17; }),
            "non-standard kernel exception did not preserve legacy fallback");
    bool contract_propagated = false;
    try {
      (void)call_legacy_occ_operation([]() -> int {
        throw stepmesher::occ::NativeContractError("[native_identity] injected");
      });
    } catch (const stepmesher::occ::NativeContractError &) {
      contract_propagated = true;
    }
    require(contract_propagated, "identity failure was hidden by numerical fallback");
    if (argc != 3) throw std::runtime_error("usage: native_contract STEP bundled-helper");
    stepmesher::MeshOptions opts;
    opts.order = 2;
    opts.mesh_fraction = .35;
    const auto nodes = stepmesher::equispaced_triangle_nodes(opts.order);
    const auto gm = stepmesher::mesh_step_file(argv[1], opts, nodes);
    require(gm.occ_evaluation == "original_gmsh", "empty selection changed numerical path");
    require(gm.native_identity_faces > 0, "identity coverage missing");
    opts.occ_info_exe = "/stepmesher-v3-test/path-that-does-not-exist";
    const auto missing = stepmesher::mesh_step_file(argv[1], opts, nodes);
    require(missing.occ_evaluation == "original_gmsh", "missing selection changed numerical path");
    require(gm.xyz == missing.xyz && gm.scaffold_triangles == missing.scaffold_triangles,
            "missing helper fallback differs from default");
    opts.occ_info_exe = argv[2];
    const auto native = stepmesher::mesh_step_file(argv[1], opts, nodes);
    require(native.occ_evaluation == "original_occ_in_process", "bundled helper not selected");
    opts.occt_precision = 1e-3;
    opts.occt_maxprecision = 1e-1;
    opts.sameparameter = false;
    opts.surfacecurve_mode = "previously_ignored_value";
    opts.cad_edge_tol_fraction = 1.4e-6; // Historical six-decimal CLI conversion.
    const auto ignored = stepmesher::mesh_step_file(argv[1], opts, nodes);
    require(native.xyz == ignored.xyz && native.panel_xyz == ignored.panel_xyz &&
            native.scaffold_triangles == ignored.scaffold_triangles &&
            native.min_meaningful_cad_edge_length == ignored.min_meaningful_cad_edge_length,
            "refactor activated old no-op settings or changed fraction rounding");
    opts.occ_info_exe = argv[0];
    bool rejected = false;
    try { (void)stepmesher::mesh_step_file(argv[1], opts, nodes); }
    catch (const std::exception &e) {
      rejected = std::string(e.what()).find("[unsupported_occ_helper]") != std::string::npos;
    }
    require(rejected, "custom executable silently accepted");
    opts.occ_info_exe = argv[2];
    const auto repeat = stepmesher::mesh_step_file(argv[1], opts, nodes);
    require(repeat.xyz == ignored.xyz, "failed selection contaminated next job");
    std::cout << "PASS native identity coverage, both numerical paths, missing/custom helpers, no-op options, rounding, repeated lifetime\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n'; return 1;
  }
}

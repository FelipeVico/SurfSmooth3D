#!/usr/bin/env python3
"""Instrument an archived mesher source copy for diagnostic benchmarks.

Never run this on the working source tree. Normal production builds contain
none of these timers or request writes. Set SS3D_PROJECTION_CAPTURE_PREFIX in
the benchmark process to collect stage times and the actual preprojection RV
points. The request/map files can be replayed by test_occ_frozen project-nodes.
"""
import argparse
from pathlib import Path


HELPER = r'''
// Test-only instrumentation injected into an archived source copy.
class ProjectionBenchmark {
  using Clock = std::chrono::steady_clock;
  std::string prefix_;
  Clock::time_point last_;
public:
  ProjectionBenchmark() : last_(Clock::now()) {
    if (const char *prefix = std::getenv("SS3D_PROJECTION_CAPTURE_PREFIX"))
      prefix_ = prefix;
    if (!prefix_.empty()) {
      std::ofstream out(prefix_ + ".timings.csv");
      out << "stage,seconds\n";
    }
  }
  void mark(const char *stage) {
    const auto now = Clock::now();
    if (!prefix_.empty()) {
      std::ofstream out(prefix_ + ".timings.csv", std::ios::app);
      out.precision(17);
      out << stage << ',' << std::chrono::duration<double>(now-last_).count() << '\n';
    }
    last_ = Clock::now();
  }
  void capture(const std::string &step_file, const occ::Context &cad,
               const std::vector<int> &face_tags,
               const std::vector<Vec3> &points, const std::vector<int> &tags) {
    if (prefix_.empty()) return;
    const auto descriptors = cad.adapter_face_descriptors();
    std::ofstream out(prefix_ + ".request.txt");
    std::ofstream map(prefix_ + ".map.txt");
    std::ofstream shape(prefix_ + ".brep");
    if (!out || !map || !shape || descriptors.size() != face_tags.size())
      throw std::runtime_error("benchmark capture initialization failed");
    out.precision(17);
    out << "STEP_MESHER_OCC_PROJECT_V1\n" << step_file << '\n' << face_tags.size() << '\n';
    for (std::size_t i = 0; i < descriptors.size(); ++i) {
      out << face_tags[i];
      for (double x : descriptors[i].bbox) out << ' ' << x;
      out << ' ' << descriptors[i].area << '\n';
      map << face_tags[i] << ' ' << i << '\n';
    }
    out << points.size() << '\n';
    for (std::size_t i = 0; i < points.size(); ++i)
      out << tags[i] << ' ' << points[i][0] << ' ' << points[i][1] << ' ' << points[i][2] << '\n';
    shape << cad.adapter_brep_snapshot();
    out.close(); map.close(); shape.close();
    if (!out || !map || !shape) throw std::runtime_error("benchmark capture write failed");
    if (std::getenv("SS3D_PROJECTION_CAPTURE_ONLY"))
      throw std::runtime_error("benchmark capture-only complete");
  }
};
'''


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError(f"Expected one instrumentation anchor: {old[:90]!r}")
    return text.replace(old, new, 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="archived SurfSmooth3D root under build/")
    args = parser.parse_args()
    source = args.source.resolve()
    live_root = Path(__file__).resolve().parents[2]
    if source == live_root or not source.is_relative_to(live_root / "build"):
        raise SystemExit("Refusing to instrument anything except a source copy under this repository's build/")
    path = source / "step_mesher/cpp/src/mesher.cpp"
    text = path.read_text()
    text = replace_once(text, "using Vec3 = std::array<double, 3>;", "using Vec3 = std::array<double, 3>;\n" + HELPER)
    anchors = [
        ("  validate_mesh_options(options);", "  validate_mesh_options(options);\n  ProjectionBenchmark benchmark;"),
        ("  auto &cad = *retained;", "  auto &cad = *retained;\n  benchmark.mark(\"import\");"),
        ("  gmsh::model::add(\"step_mesher_cpp_matlab\");", "  benchmark.mark(\"cad_measurements\");\n  gmsh::model::add(\"step_mesher_cpp_matlab\");"),
        ("  const auto fragmentation = apply_continuity_fragmentation(", "  benchmark.mark(\"binding_and_bbox\");\n  const auto fragmentation = apply_continuity_fragmentation("),
        ("  configure_meshing(result.mesh_size, result.hmin_effective,", "  benchmark.mark(\"fragmentation\");\n  configure_meshing(result.mesh_size, result.hmin_effective,"),
        ("  gmsh::model::mesh::generate(2);", "  gmsh::model::mesh::generate(2);\n  benchmark.mark(\"scaffold\");"),
        ("  std::vector<Vec3> candidate_xyz;", "  benchmark.mark(\"edge_processing\");\n  std::vector<Vec3> candidate_xyz;"),
        ("  auto projected = load_occ_projected_nodes(", "  benchmark.mark(\"candidate_generation\");\n  benchmark.capture(step_file, cad, bound_face_tags, candidate_xyz, candidate_face_tags);\n  benchmark.mark(\"diagnostic_capture\");\n  auto projected = load_occ_projected_nodes("),
        ("  if (projected.ok) {", "  benchmark.mark(\"projection\");\n  if (projected.ok) {"),
        ("      projection_count > 0 ? sum_projection_distance / projection_count : 0.0;", "      projection_count > 0 ? sum_projection_distance / projection_count : 0.0;\n  benchmark.mark(\"result_packaging\");"),
    ]
    for old, new in anchors:
        text = replace_once(text, old, new)
    path.write_text(text)
    print(path)


if __name__ == "__main__":
    main()

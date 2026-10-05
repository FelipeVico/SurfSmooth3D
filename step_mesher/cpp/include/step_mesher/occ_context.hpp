#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <string>

namespace stepmesher::occ {

// Internal integration contract failures must not become legacy numerical
// fallbacks at the in-process replacement for the helper protocol boundary.
class NativeContractError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

using Vec3 = std::array<double, 3>;
using Vec2 = std::array<double, 2>;

struct Line { std::size_t a = 0, b = 0; int curve = 0; };
struct EdgeRequest {
  std::vector<double> edge_nodes;
  std::unordered_map<std::size_t, Vec3> nodes;
  std::unordered_map<int, std::vector<int>> curve_faces;
  std::vector<Line> lines;
};
struct LineSamples { std::size_t a = 0, b = 0; std::vector<Vec3> points; };
struct EdgeSamples {
  std::vector<LineSamples> lines;
  int n_newton = 0, n_curve = 0, n_fallback = 0;
};
struct HelperCurve {
  int face_tag = 0, axis = -1;
  double param = std::numeric_limits<double>::quiet_NaN();
  std::vector<Vec3> points;
};
struct HelperCurves { std::vector<HelperCurve> curves; int break_count = 0; };
struct FeatureLine {
  std::size_t a = 0, b = 0;
  int face_tag = 0, axis = -1;
  double param = std::numeric_limits<double>::quiet_NaN();
};
struct FeatureRequest {
  std::vector<double> edge_nodes;
  std::unordered_map<std::size_t, Vec3> nodes;
  std::vector<FeatureLine> lines;
};
struct FeatureSamples { std::vector<LineSamples> lines; int n_fallback = 0; };
struct FeatureNodeConstraint {
  int face_tag = 0, axis = -1;
  double param = std::numeric_limits<double>::quiet_NaN();
};
struct FeatureNode {
  std::size_t tag = 0;
  Vec3 x{};
  std::vector<FeatureNodeConstraint> constraints;
};
struct CorrectedFeatureNode { std::size_t tag = 0; Vec3 x{}; double distance = 0.0; };
struct FeatureNodes { std::vector<CorrectedFeatureNode> nodes; int n_fallback = 0; };
struct ProjectNode { int face_tag = 0; Vec3 x{}; };
struct ProjectedNode { Vec3 x{}; double distance = 0.0; Vec2 uv{}; };
struct ProjectedNodes { std::vector<ProjectedNode> nodes; int n_fallback = 0; };
struct Measurements {
  double diag = 0.0, area = 0.0;
  double min_raw_edge = std::numeric_limits<double>::quiet_NaN();
  double min_meaningful_edge = std::numeric_limits<double>::quiet_NaN();
  std::size_t n_edges = 0;
  int n_meaningful = 0;
};

// Descriptors exist solely for the standalone protocol adapter and differential
// tests. The production native bridge never compares these quantities.
struct FaceDescriptor {
  int tag = 0;
  std::array<double, 6> bbox{};
  Vec3 center{}, size{};
  double diag = 0.0, area = 0.0;
};
struct FaceOccurrence { std::size_t canonical_index = 0; int orientation = 0; };

// A context owns one retained STEP import. Gmsh may augment edge p-curve
// representations during fragmentation. Its lifetime encloses all Gmsh
// operations on that imported shape. Caller serializes process-global CAD use.
class Context {
public:
  explicit Context(const std::string &step_file);
  ~Context();
  Context(const Context &) = delete;
  Context &operator=(const Context &) = delete;
  Context(Context &&) noexcept;
  Context &operator=(Context &&) noexcept;

  // Imports the complete retained shape, then resolves its retained faces by
  // native topology identity. Gmsh must already have a current empty model.
  std::vector<int> bind_native();
  std::size_t face_count() const;
  const std::vector<FaceOccurrence> &face_occurrences() const;
  Measurements measurements(double cad_edge_tol_fraction) const;
  HelperCurves helper_curves(double target_spacing,
                            const std::vector<int> &ordered_source_tags) const;
  EdgeSamples edge_samples(const EdgeRequest &req) const;
  FeatureNodes feature_nodes(const std::vector<FeatureNode> &nodes) const;
  FeatureSamples feature_samples(const FeatureRequest &req) const;
  ProjectedNodes project_nodes(const std::vector<ProjectNode> &nodes) const;

  // Adapter/test-only binding: zero-based canonical face indices. Not called
  // by native production construction; exact bindings otherwise come above.
  std::vector<FaceDescriptor> adapter_face_descriptors() const;
  // Test-only complete BRep evidence, omitting triangulations.
  std::string adapter_brep_snapshot() const;
  void bind_adapter_face_indices(const std::unordered_map<int, int> &indices);

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace stepmesher::occ

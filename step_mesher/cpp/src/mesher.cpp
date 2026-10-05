#include "step_mesher/mesher.hpp"
#include "step_mesher/occ_context.hpp"
#include "occ_operation.hpp"
#include <memory>
#include <mutex>
#include <dlfcn.h>
#include <Standard.hxx>
#include <Standard_Version.hxx>

#include <gmsh.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace stepmesher {
namespace {

using Vec2 = std::array<double, 2>;
using Vec3 = std::array<double, 3>;

struct NodeRecord {
  Vec3 xyz;
  int local_index = 0;
};

struct EdgeKey {
  std::size_t a = 0;
  std::size_t b = 0;

  bool operator==(const EdgeKey &other) const {
    return a == other.a && b == other.b;
  }
};

struct EdgeKeyHash {
  std::size_t operator()(const EdgeKey &key) const {
    return std::hash<std::size_t>{}(key.a) ^
           (std::hash<std::size_t>{}(key.b) + 0x9e3779b97f4a7c15ULL +
            (std::hash<std::size_t>{}(key.a) << 6) +
            (std::hash<std::size_t>{}(key.a) >> 2));
  }
};

using EdgeSampleMap =
    std::unordered_map<EdgeKey, std::vector<Vec3>, EdgeKeyHash>;
using CurveFaceMap = std::unordered_map<int, std::vector<int>>;

struct OccInfo {
  bool ok = false;
  double tight_bbox_diag = std::numeric_limits<double>::quiet_NaN();
  double surface_area = std::numeric_limits<double>::quiet_NaN();
  double min_raw_cad_edge_length = std::numeric_limits<double>::quiet_NaN();
  double min_meaningful_cad_edge_length =
      std::numeric_limits<double>::quiet_NaN();
  int cad_edge_count = 0;
  int meaningful_cad_edge_count = 0;
};

struct OccEdgeHelperResult {
  bool ok = false;
  EdgeSampleMap edge_samples;
  int newton_success = 0;
  int curve_projection_success = 0;
  int fallback_points = 0;
};

struct OccProjectHelperResult {
  bool ok = false;
  std::vector<Vec3> xyz;
  std::vector<double> dist;
  int fallback_points = 0;
};

struct OccFeatureNodeHelperResult {
  bool ok = false;
  std::unordered_map<std::size_t, Vec3> nodes;
  int fallback_points = 0;
};

struct TriangleQuality {
  double min_angle_deg = 0.0;
  double edge_ratio = std::numeric_limits<double>::infinity();
  double area = 0.0;
};

struct ContinuityHelperCurve {
  int original_face_tag = 0;
  int axis = -1;
  double param = std::numeric_limits<double>::quiet_NaN();
  std::vector<Vec3> points;
};

struct OccContinuityHelperResult {
  bool ok = false;
  std::vector<ContinuityHelperCurve> curves;
  int break_count = 0;
};

struct FeatureCurveInfo {
  int original_face_tag = 0;
  int axis = -1;
  double param = std::numeric_limits<double>::quiet_NaN();
};

using FeatureCurveMap = std::unordered_map<int, FeatureCurveInfo>;

struct FeatureNodeConstraint {
  int original_face_tag = 0;
  int axis = -1;
  double param = std::numeric_limits<double>::quiet_NaN();
};

using FeatureNodeConstraintMap =
    std::unordered_map<std::size_t, std::vector<FeatureNodeConstraint>>;

struct FragmentationInfo {
  std::unordered_map<int, int> current_face_to_original_face;
  FeatureCurveMap feature_curves;
  int helper_curve_count = 0;
  int continuity_break_count = 0;
  int fragmented_face_count = 0;
};

EdgeKey unoriented_edge(std::size_t a, std::size_t b) {
  return a < b ? EdgeKey{a, b} : EdgeKey{b, a};
}

void append_unique(std::vector<int> &values, int value) {
  if (std::find(values.begin(), values.end(), value) == values.end()) {
    values.push_back(value);
  }
}

class GmshSession {
public:
  explicit GmshSession(bool verbose) {
    was_initialized_ = gmsh::isInitialized() != 0;
    if (!was_initialized_) {
      gmsh::initialize();
    }
    gmsh::option::setNumber("General.Terminal", verbose ? 1.0 : 0.0);
    gmsh::clear();
  }

  ~GmshSession() {
    try {
      gmsh::clear();
      if (!was_initialized_) {
        gmsh::finalize();
      }
    } catch (...) {
    }
  }

  GmshSession(const GmshSession &) = delete;
  GmshSession &operator=(const GmshSession &) = delete;

private:
  bool was_initialized_ = false;
};

void set_option_if_available(const std::string &name, double value) {
  try {
    gmsh::option::setNumber(name, value);
  } catch (...) {
  }
}

double norm(const Vec3 &a) {
  return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

double clamp_unit(double x) {
  return std::max(-1.0, std::min(1.0, x));
}

bool is_finite(const Vec3 &a) {
  return std::isfinite(a[0]) && std::isfinite(a[1]) && std::isfinite(a[2]);
}

double dot(const Vec3 &a, const Vec3 &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 cross(const Vec3 &a, const Vec3 &b) {
  return {a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

Vec3 subtract(const Vec3 &a, const Vec3 &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec3 blend3_bary(const Vec3 &a, const Vec3 &b, const Vec3 &c,
                 const Vec3 &bary) {
  return {bary[0] * a[0] + bary[1] * b[0] + bary[2] * c[0],
          bary[0] * a[1] + bary[1] * b[1] + bary[2] * c[1],
          bary[0] * a[2] + bary[1] * b[2] + bary[2] * c[2]};
}

Vec3 blend_segment(const Vec3 &a, const Vec3 &b, double t) {
  return {(1.0 - t) * a[0] + t * b[0],
          (1.0 - t) * a[1] + t * b[1],
          (1.0 - t) * a[2] + t * b[2]};
}

Vec3 add(const Vec3 &a, const Vec3 &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

TriangleQuality compute_triangle_quality(const Vec3 &p0, const Vec3 &p1,
                                         const Vec3 &p2) {
  TriangleQuality q;
  const Vec3 e01 = subtract(p1, p0);
  const Vec3 e12 = subtract(p2, p1);
  const Vec3 e20 = subtract(p0, p2);
  const double l01 = norm(e01);
  const double l12 = norm(e12);
  const double l20 = norm(e20);
  const double lmin = std::min({l01, l12, l20});
  const double lmax = std::max({l01, l12, l20});
  q.edge_ratio = lmin > 0.0 ? lmax / lmin
                            : std::numeric_limits<double>::infinity();
  q.area = 0.5 * norm(cross(e01, subtract(p2, p0)));
  if (l01 > 0.0 && l12 > 0.0 && l20 > 0.0) {
    constexpr double rad_to_deg = 57.2957795130823208768;
    const double a0 =
        std::acos(clamp_unit(dot(e01, subtract(p2, p0)) / (l01 * l20))) *
        rad_to_deg;
    const double a1 =
        std::acos(clamp_unit(dot(subtract(p0, p1), e12) / (l01 * l12))) *
        rad_to_deg;
    const double a2 =
        std::acos(clamp_unit(dot(subtract(p0, p2), subtract(p1, p2)) /
                             (l20 * l12))) *
        rad_to_deg;
    q.min_angle_deg = std::min({a0, a1, a2});
  }
  return q;
}

Vec3 scale(const Vec3 &a, double s) {
  return {s * a[0], s * a[1], s * a[2]};
}

// The six former subprocess calls now use the retained native CAD context.
// Numerical choices and the historical enabled/disabled split are preserved.
std::mutex cad_process_mutex;

void verify_native_runtime() {
  const std::filesystem::path expected =
      std::filesystem::canonical(STEP_MESHER_DEPENDENCY_PREFIX) / "lib";
  const std::pair<const char *, const void *> libraries[] = {
    {"Gmsh", reinterpret_cast<const void *>(&gmsh::initialize)},
    {"OCCT", reinterpret_cast<const void *>(&Standard::Allocate)}
  };
  for (const auto &library : libraries) {
    Dl_info loaded{};
    if (!dladdr(library.second, &loaded) || !loaded.dli_fname ||
        std::filesystem::canonical(loaded.dli_fname).parent_path() != expected)
      throw std::runtime_error(std::string("[native_runtime] ") + library.first +
          " did not load from V3's common dependency installation");
  }
  std::string version;
  gmsh::option::getString("General.Version", version);
  if (version != "4.12.2" || std::string(OCC_VERSION_COMPLETE) != "7.9.3")
    throw std::runtime_error("[native_runtime] Expected Gmsh 4.12.2 / OCCT 7.9.3");
}

bool same_file_contents(const std::filesystem::path &a,
                        const std::filesystem::path &b) {
  std::error_code ec;
  if (!std::filesystem::is_regular_file(b, ec) ||
      std::filesystem::file_size(a, ec) != std::filesystem::file_size(b, ec))
    return false;
  if (std::filesystem::equivalent(a, b, ec)) return true;
  std::ifstream x(a, std::ios::binary), y(b, std::ios::binary);
  char xb[16384], yb[16384];
  while (x && y) {
    x.read(xb, sizeof(xb)); y.read(yb, sizeof(yb));
    if (x.gcount() != y.gcount() ||
        !std::equal(xb, xb + x.gcount(), yb)) return false;
  }
  return x.eof() && y.eof();
}

bool use_occ_computations(const std::string &helper) {
  if (helper.empty()) return false;
  std::error_code ec;
  if (!std::filesystem::exists(helper, ec)) return false;
  Dl_info loaded{};
  if (!dladdr(reinterpret_cast<void *>(&use_occ_computations), &loaded) ||
      !loaded.dli_fname)
    throw std::runtime_error("[native_runtime] Cannot locate the loaded v3 core");
  auto libdir = std::filesystem::canonical(loaded.dli_fname).parent_path();
  const auto bundled = libdir.parent_path() /
      "matlab/+surfsmooth3d/+stepmesher/private/step_mesher_occ_info";
  if (!same_file_contents(helper, bundled))
    throw std::runtime_error("[unsupported_occ_helper] Use this v3 runtime's bundled step_mesher_occ_info; custom helpers are not executed: " + helper);
  return true;
}

OccInfo load_occ_info(occ::Context &cad, bool enabled,
                      double cad_edge_tol_fraction) {
  OccInfo info;
  if (!enabled) return info;
  // The external helper received this argument with six decimal places.
  const double compatibility_fraction = std::stod(std::to_string(cad_edge_tol_fraction));
  const auto operation = detail::call_legacy_occ_operation(
      [&] { return cad.measurements(compatibility_fraction); });
  if (!operation) return info;
  const auto &measured = *operation;
  info.tight_bbox_diag = measured.diag;
  info.surface_area = measured.area;
  info.min_raw_cad_edge_length = measured.min_raw_edge;
  info.min_meaningful_cad_edge_length = measured.min_meaningful_edge;
  info.cad_edge_count = static_cast<int>(measured.n_edges);
  info.meaningful_cad_edge_count = measured.n_meaningful;
  info.ok = std::isfinite(info.tight_bbox_diag) && info.tight_bbox_diag > 0.0 &&
            std::isfinite(info.surface_area) && info.surface_area > 0.0;
  return info;
}

bool store_line_samples(const std::vector<occ::LineSamples> &lines,
                        std::size_t count, EdgeSampleMap &out) {
  for (const auto &line : lines) {
    if (line.points.size() != count ||
        !std::all_of(line.points.begin(), line.points.end(), is_finite)) return false;
    out[EdgeKey{line.a,line.b}] = line.points;
    std::vector<Vec3> reversed = line.points;
    std::reverse(reversed.begin(), reversed.end());
    out[EdgeKey{line.b,line.a}] = std::move(reversed);
  }
  return true;
}

OccEdgeHelperResult load_occ_edge_samples(
    occ::Context &cad, bool enabled,
    const std::vector<double> &edge_nodes, const CurveFaceMap &curve_faces,
    const std::unordered_map<std::size_t, NodeRecord> &node_map,
    const FeatureCurveMap &feature_curves) {
  OccEdgeHelperResult result;
  if (!enabled) return result;
  occ::EdgeRequest request;
  request.edge_nodes = edge_nodes;
  for (const auto &entry : node_map) request.nodes[entry.first] = entry.second.xyz;
  for (const auto &entry : curve_faces)
    if (feature_curves.find(entry.first) == feature_curves.end())
      request.curve_faces[entry.first] = entry.second;
  std::vector<occ::Line> lines;
  gmsh::vectorpair curves;
  gmsh::model::getEntities(curves, 1);
  for (const auto &curve : curves) {
    const int curve_tag = curve.second;
    if (feature_curves.find(curve_tag) != feature_curves.end()) {
      continue;
    }
    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;
    gmsh::model::mesh::getElements(element_types, element_tags,
                                   element_node_tags, 1, curve_tag);
    for (std::size_t block = 0; block < element_types.size(); ++block) {
      std::string element_name;
      int dim = 0;
      int element_order = 0;
      int num_nodes = 0;
      int num_primary_nodes = 0;
      std::vector<double> local_coords;
      gmsh::model::mesh::getElementProperties(
          element_types[block], element_name, dim, element_order, num_nodes,
          local_coords, num_primary_nodes);
      if (dim != 1 || num_primary_nodes < 2 ||
          element_name.find("Line") == std::string::npos) {
        continue;
      }
      const auto &tags = element_node_tags[block];
      const std::size_t num_elements =
          tags.size() / static_cast<std::size_t>(num_nodes);
      for (std::size_t elem = 0; elem < num_elements; ++elem) {
        lines.push_back({tags[elem * static_cast<std::size_t>(num_nodes)],
                         tags[elem * static_cast<std::size_t>(num_nodes) + 1],
                         curve_tag});
      }
    }
  }

  request.lines = std::move(lines);
  const auto operation = detail::call_legacy_occ_operation(
      [&] { return cad.edge_samples(request); });
  if (!operation) return result;
  const auto &sampled = *operation;
  if (sampled.lines.empty()) return result;
  result.newton_success = sampled.n_newton;
  result.curve_projection_success = sampled.n_curve;
  result.fallback_points = sampled.n_fallback;
  if (!store_line_samples(sampled.lines, edge_nodes.size(), result.edge_samples))
    return OccEdgeHelperResult{};
  result.ok = true;
  return result;
}

OccContinuityHelperResult load_occ_helper_curves(
    occ::Context &cad, bool enabled, double target_spacing,
    const std::vector<int> &source_tags) {
  OccContinuityHelperResult result;
  if (!enabled || source_tags.empty()) return result;
  const auto operation = detail::call_legacy_occ_operation(
      [&] { return cad.helper_curves(target_spacing, source_tags); });
  if (!operation) return result;
  const auto &helpers = *operation;
  result.break_count = helpers.break_count;
  for (const auto &curve : helpers.curves) {
    if (curve.points.size() < 2 ||
        !std::all_of(curve.points.begin(), curve.points.end(), is_finite))
      return OccContinuityHelperResult{};
    result.curves.push_back({curve.face_tag, curve.axis, curve.param, curve.points});
  }
  result.ok = true;
  return result;
}

OccEdgeHelperResult load_occ_feature_samples(
    occ::Context &cad, bool enabled, const std::vector<double> &edge_nodes,
    const std::unordered_map<std::size_t, NodeRecord> &node_map,
    const FeatureCurveMap &feature_curves) {
  OccEdgeHelperResult result;
  if (!enabled || feature_curves.empty()) return result;
  occ::FeatureRequest request;
  request.edge_nodes = edge_nodes;
  for (const auto &entry : node_map) request.nodes[entry.first] = entry.second.xyz;
  struct FeatureSegment {
    std::size_t a = 0;
    std::size_t b = 0;
    FeatureCurveInfo info;
  };
  std::vector<FeatureSegment> segments;
  for (const auto &entry : feature_curves) {
    const int curve_tag = entry.first;
    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;
    gmsh::model::mesh::getElements(element_types, element_tags,
                                   element_node_tags, 1, curve_tag);
    for (std::size_t block = 0; block < element_types.size(); ++block) {
      std::string element_name;
      int dim = 0;
      int element_order = 0;
      int num_nodes = 0;
      int num_primary_nodes = 0;
      std::vector<double> local_coords;
      gmsh::model::mesh::getElementProperties(
          element_types[block], element_name, dim, element_order, num_nodes,
          local_coords, num_primary_nodes);
      if (dim != 1 || num_primary_nodes < 2 ||
          element_name.find("Line") == std::string::npos) {
        continue;
      }
      const auto &tags = element_node_tags[block];
      const std::size_t num_elements =
          tags.size() / static_cast<std::size_t>(num_nodes);
      for (std::size_t elem = 0; elem < num_elements; ++elem) {
        segments.push_back({tags[elem * static_cast<std::size_t>(num_nodes)],
                            tags[elem * static_cast<std::size_t>(num_nodes) + 1],
                            entry.second});
      }
    }
  }

  for (const auto &segment : segments)
    request.lines.push_back({segment.a, segment.b, segment.info.original_face_tag,
                             segment.info.axis, segment.info.param});
  const auto operation = detail::call_legacy_occ_operation(
      [&] { return cad.feature_samples(request); });
  if (!operation) return result;
  const auto &sampled = *operation;
  result.fallback_points = sampled.n_fallback;
  if (!store_line_samples(sampled.lines, edge_nodes.size(), result.edge_samples))
    return OccEdgeHelperResult{};
  result.ok = true;
  return result;
}

OccFeatureNodeHelperResult load_occ_feature_nodes(
    occ::Context &cad, bool enabled,
    const std::unordered_map<std::size_t, Vec3> &raw_nodes,
    const FeatureNodeConstraintMap &constraints) {
  OccFeatureNodeHelperResult result;
  if (!enabled || constraints.empty()) return result;
  std::vector<occ::FeatureNode> request;
  for (const auto &entry : constraints) {
    const auto it = raw_nodes.find(entry.first);
    if (it == raw_nodes.end() || entry.second.empty()) continue;
    occ::FeatureNode node;
    node.tag = entry.first; node.x = it->second;
    for (const auto &constraint : entry.second)
      node.constraints.push_back({constraint.original_face_tag, constraint.axis,
                                 constraint.param});
    request.push_back(std::move(node));
  }
  const auto operation = detail::call_legacy_occ_operation(
      [&] { return cad.feature_nodes(request); });
  if (!operation) return result;
  const auto &corrected = *operation;
  result.fallback_points = corrected.n_fallback;
  for (const auto &node : corrected.nodes) {
    if (!is_finite(node.x) || !std::isfinite(node.distance))
      return OccFeatureNodeHelperResult{};
    result.nodes[node.tag] = node.x;
  }
  result.ok = true;
  return result;
}

OccProjectHelperResult load_occ_projected_nodes(
    occ::Context &cad, bool enabled, const std::vector<Vec3> &candidate_xyz,
    const std::vector<int> &candidate_face_tags) {
  OccProjectHelperResult result;
  if (candidate_xyz.empty()) { result.ok = true; return result; }
  if (!enabled) return result;
  if (candidate_xyz.size() != candidate_face_tags.size()) return result;
  std::vector<occ::ProjectNode> request;
  request.reserve(candidate_xyz.size());
  for (std::size_t i = 0; i < candidate_xyz.size(); ++i)
    request.push_back({candidate_face_tags[i], candidate_xyz[i]});
  const auto operation = detail::call_legacy_occ_operation(
      [&] { return cad.project_nodes(request); });
  if (!operation) return result;
  const auto &projected = *operation;
  result.fallback_points = projected.n_fallback;
  if (projected.nodes.size() != candidate_xyz.size()) return OccProjectHelperResult{};
  result.xyz.reserve(projected.nodes.size()); result.dist.reserve(projected.nodes.size());
  for (const auto &node : projected.nodes) {
    if (!is_finite(node.x) || !std::isfinite(node.distance)) return OccProjectHelperResult{};
    result.xyz.push_back(node.x); result.dist.push_back(node.distance);
  }
  result.ok = true;
  return result;
}

struct RefTriangle {
  Vec3 a;
  Vec3 b;
  Vec3 c;
};

std::vector<RefTriangle> refined_reference_triangles(int refinement_level) {
  if (refinement_level < 0) {
    throw std::invalid_argument("refinement_level must be nonnegative");
  }

  const int n = 1 << refinement_level;
  std::vector<RefTriangle> tris;
  tris.reserve(static_cast<std::size_t>(n * n));

  const auto bary = [n](int i, int j) -> Vec3 {
    const double u = static_cast<double>(i) / static_cast<double>(n);
    const double v = static_cast<double>(j) / static_cast<double>(n);
    // Store true barycentric weights relative to (P1, P2, P3). The
    // reference coordinates (u,v) follow the fmm3dbie/koorn convention:
    // x(u,v) = (1-u-v) P1 + u P2 + v P3.
    return {1.0 - u - v, u, v};
  };

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n - i; ++j) {
      const Vec3 p00 = bary(i, j);
      const Vec3 p10 = bary(i + 1, j);
      const Vec3 p01 = bary(i, j + 1);
      tris.push_back({p00, p10, p01});

      if (i + j < n - 1) {
        const Vec3 p11 = bary(i + 1, j + 1);
        tris.push_back({p10, p11, p01});
      }
    }
  }
  return tris;
}

Vec3 map_local_bary_to_parent(const Vec3 &local_bary,
                              const RefTriangle &subtri) {
  return {local_bary[0] * subtri.a[0] + local_bary[1] * subtri.b[0] +
              local_bary[2] * subtri.c[0],
          local_bary[0] * subtri.a[1] + local_bary[1] * subtri.b[1] +
              local_bary[2] * subtri.c[1],
          local_bary[0] * subtri.a[2] + local_bary[1] * subtri.b[2] +
              local_bary[2] * subtri.c[2]};
}

void legendre_with_derivative(int n, double x, double &pn, double &dpn) {
  if (n == 0) {
    pn = 1.0;
    dpn = 0.0;
    return;
  }
  double pnm2 = 1.0;
  double pnm1 = x;
  if (n == 1) {
    pn = pnm1;
  } else {
    for (int k = 2; k <= n; ++k) {
      const double pk =
          ((2.0 * k - 1.0) * x * pnm1 - (k - 1.0) * pnm2) /
          static_cast<double>(k);
      pnm2 = pnm1;
      pnm1 = pk;
    }
    pn = pnm1;
  }

  const double den = x * x - 1.0;
  if (std::abs(den) < 1.0e-14) {
    dpn = 0.5 * n * (n + 1) * (x >= 0.0 ? 1.0 : std::pow(-1.0, n - 1));
  } else {
    dpn = static_cast<double>(n) * (x * pn - pnm2) / den;
  }
}

std::vector<double> legendre_gauss_lobatto_nodes(int order) {
  const int n = std::max(order, 1);
  std::vector<double> nodes(static_cast<std::size_t>(n + 1));
  nodes.front() = 0.0;
  nodes.back() = 1.0;
  if (n == 1) {
    return nodes;
  }

  constexpr double pi = 3.141592653589793238462643383279502884;
  for (int i = 1; i < n; ++i) {
    // Interior Lobatto nodes are the roots of P'_n. The cosine seed is the
    // standard asymptotic starting point; mapping to [0, 1] happens below.
    double x = -std::cos(pi * static_cast<double>(i) / static_cast<double>(n));
    for (int iter = 0; iter < 50; ++iter) {
      double pn = 0.0;
      double dpn = 0.0;
      legendre_with_derivative(n, x, pn, dpn);

      // Legendre's ODE gives P''_n at a root solve point.
      const double denom = 1.0 - x * x;
      if (std::abs(denom) < 1.0e-15) {
        break;
      }
      const double ddpn =
          (2.0 * x * dpn - n * (n + 1.0) * pn) / denom;
      const double dx = dpn / ddpn;
      x -= dx;
      if (std::abs(dx) < 4.0e-16) {
        break;
      }
    }
    nodes[static_cast<std::size_t>(i)] = 0.5 * (x + 1.0);
  }
  return nodes;
}

std::vector<double>
legendre_gauss_lobatto_quadrature_weights(const std::vector<double> &nodes) {
  const int n = static_cast<int>(nodes.size()) - 1;
  if (n <= 0) {
    return {1.0};
  }

  std::vector<double> weights(nodes.size(), 0.0);
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const double x = 2.0 * nodes[i] - 1.0;
    double pn = 0.0;
    double dpn = 0.0;
    legendre_with_derivative(n, x, pn, dpn);
    weights[i] = 1.0 / (static_cast<double>(n) * (n + 1.0) * pn * pn);
  }
  return weights;
}

std::vector<double> barycentric_weights(const std::vector<double> &nodes) {
  std::vector<double> weights(nodes.size(), 1.0);
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    for (std::size_t k = 0; k < nodes.size(); ++k) {
      if (j != k) {
        weights[j] /= nodes[j] - nodes[k];
      }
    }
  }
  return weights;
}

Vec3 interpolate_edge(const std::vector<double> &nodes,
                      const std::vector<double> &weights,
                      const std::vector<Vec3> &points, double t) {
  constexpr double tol = 1.0e-14;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (std::abs(t - nodes[i]) <= tol) {
      return points[i];
    }
  }

  Vec3 numerator{0.0, 0.0, 0.0};
  double denominator = 0.0;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const double tmp = weights[i] / (t - nodes[i]);
    numerator = add(numerator, scale(points[i], tmp));
    denominator += tmp;
  }
  return scale(numerator, 1.0 / denominator);
}

std::vector<Vec3> edge_derivatives_at_nodes(
    const std::vector<double> &nodes, const std::vector<double> &weights,
    const std::vector<Vec3> &points) {
  std::vector<Vec3> derivs(points.size(), Vec3{0.0, 0.0, 0.0});
  if (nodes.size() != points.size() || weights.size() != points.size()) {
    return derivs;
  }

  for (std::size_t i = 0; i < points.size(); ++i) {
    Vec3 d{0.0, 0.0, 0.0};
    for (std::size_t j = 0; j < points.size(); ++j) {
      if (i == j) {
        continue;
      }
      const double denom = nodes[i] - nodes[j];
      if (std::abs(denom) < 1.0e-300 || weights[i] == 0.0) {
        continue;
      }
      const double coeff = weights[j] / (weights[i] * denom);
      d = add(d, scale(subtract(points[j], points[i]), coeff));
    }
    derivs[i] = d;
  }
  return derivs;
}

Vec3 finite_difference_edge_tangent(const std::vector<double> &nodes,
                                    const std::vector<Vec3> &points,
                                    std::size_t i) {
  if (points.size() < 2 || nodes.size() != points.size()) {
    return {1.0, 0.0, 0.0};
  }

  std::size_t lo = 0;
  std::size_t hi = 0;
  if (i == 0) {
    lo = 0;
    hi = 1;
  } else if (i + 1 == points.size()) {
    lo = points.size() - 2;
    hi = points.size() - 1;
  } else {
    lo = i - 1;
    hi = i + 1;
  }

  const double da = nodes[hi] - nodes[lo];
  if (std::abs(da) < 1.0e-300) {
    return subtract(points.back(), points.front());
  }
  return scale(subtract(points[hi], points[lo]), 1.0 / da);
}

double snap01(double x) {
  constexpr double tol = 1.0e-14;
  if (std::abs(x) < tol) {
    return 0.0;
  }
  if (std::abs(x - 1.0) < tol) {
    return 1.0;
  }
  return x;
}

Vec3 evaluate_blended_triangle_point(
    Vec3 bary, const Vec3 &v1, const Vec3 &v2, const Vec3 &v3,
    const std::vector<double> &edge_nodes,
    const std::vector<double> &edge_weights,
    const std::vector<Vec3> *edge12, const std::vector<Vec3> *edge23,
    const std::vector<Vec3> *edge31) {
  constexpr double tol = 1.0e-14;
  bary[0] = snap01(bary[0]);
  bary[1] = snap01(bary[1]);
  bary[2] = snap01(bary[2]);

  const double u = bary[0];
  const double v = bary[1];
  const double w = bary[2];

  if (std::abs(u - 1.0) < tol) {
    return v1;
  }
  if (std::abs(v - 1.0) < tol) {
    return v2;
  }
  if (std::abs(w - 1.0) < tol) {
    return v3;
  }

  if (std::abs(w) < tol) {
    return edge12 ? interpolate_edge(edge_nodes, edge_weights, *edge12, v)
                  : blend_segment(v1, v2, v);
  }
  if (std::abs(u) < tol) {
    return edge23 ? interpolate_edge(edge_nodes, edge_weights, *edge23, w)
                  : blend_segment(v2, v3, w);
  }
  if (std::abs(v) < tol) {
    return edge31 ? interpolate_edge(edge_nodes, edge_weights, *edge31, u)
                  : blend_segment(v3, v1, u);
  }

  Vec3 point = blend3_bary(v1, v2, v3, bary);

  if (edge12) {
    const Vec3 curve = interpolate_edge(edge_nodes, edge_weights, *edge12, v);
    const Vec3 linear = blend_segment(v1, v2, v);
    const double factor = (u * v) / std::max(v * (1.0 - v), 1.0e-300);
    point = add(point, scale(subtract(curve, linear), factor));
  }
  if (edge23) {
    const Vec3 curve = interpolate_edge(edge_nodes, edge_weights, *edge23, w);
    const Vec3 linear = blend_segment(v2, v3, w);
    const double factor = (v * w) / std::max(w * (1.0 - w), 1.0e-300);
    point = add(point, scale(subtract(curve, linear), factor));
  }
  if (edge31) {
    const Vec3 curve = interpolate_edge(edge_nodes, edge_weights, *edge31, u);
    const Vec3 linear = blend_segment(v3, v1, u);
    const double factor = (w * u) / std::max(u * (1.0 - u), 1.0e-300);
    point = add(point, scale(subtract(curve, linear), factor));
  }

  return point;
}

bool surface_uv_from_xyz(int face_tag, const Vec3 &xyz, Vec2 &uv) {
  try {
    std::vector<double> coords{xyz[0], xyz[1], xyz[2]};
    std::vector<double> param;
    gmsh::model::getParametrization(2, face_tag, coords, param);
    if (param.size() < 2 || !std::isfinite(param[0]) ||
        !std::isfinite(param[1])) {
      return false;
    }
    uv = {param[0], param[1]};
    return true;
  } catch (...) {
    return false;
  }
}

bool xyz_from_surface_uv(int face_tag, const Vec2 &uv, Vec3 &xyz) {
  try {
    std::vector<double> param{uv[0], uv[1]};
    std::vector<double> coords;
    gmsh::model::getValue(2, face_tag, param, coords);
    if (coords.size() < 3 || !std::isfinite(coords[0]) ||
        !std::isfinite(coords[1]) || !std::isfinite(coords[2])) {
      return false;
    }
    xyz = {coords[0], coords[1], coords[2]};
    return true;
  } catch (...) {
    return false;
  }
}

bool closest_point_on_surface(int face_tag, const Vec3 &xyz, Vec3 &projected,
                              Vec2 &uv) {
  try {
    std::vector<double> coords{xyz[0], xyz[1], xyz[2]};
    std::vector<double> closest;
    std::vector<double> param;
    gmsh::model::getClosestPoint(2, face_tag, coords, closest, param);
    if (closest.size() < 3 || param.size() < 2 ||
        !std::isfinite(closest[0]) || !std::isfinite(closest[1]) ||
        !std::isfinite(closest[2]) || !std::isfinite(param[0]) ||
        !std::isfinite(param[1])) {
      return false;
    }
    projected = {closest[0], closest[1], closest[2]};
    uv = {param[0], param[1]};
    return true;
  } catch (...) {
    return false;
  }
}

bool closest_point_on_curve(int curve_tag, const Vec3 &xyz, Vec3 &projected,
                            double *curve_param = nullptr) {
  try {
    std::vector<double> coords{xyz[0], xyz[1], xyz[2]};
    std::vector<double> closest;
    std::vector<double> param;
    gmsh::model::getClosestPoint(1, curve_tag, coords, closest, param);
    if (closest.size() < 3 || !std::isfinite(closest[0]) ||
        !std::isfinite(closest[1]) || !std::isfinite(closest[2])) {
      return false;
    }
    projected = {closest[0], closest[1], closest[2]};
    if (curve_param != nullptr && !param.empty() && std::isfinite(param[0])) {
      *curve_param = param[0];
    }
    return true;
  } catch (...) {
    return false;
  }
}

struct SurfaceEval {
  Vec3 x{};
  Vec3 xu{};
  Vec3 xv{};
  Vec3 xuu{};
  Vec3 xvv{};
  Vec3 xuv{};
};

bool surface_eval_d2(int face_tag, const Vec2 &uv, SurfaceEval &eval) {
  try {
    const std::vector<double> param{uv[0], uv[1]};
    std::vector<double> coords;
    std::vector<double> d1;
    std::vector<double> d2;
    gmsh::model::getValue(2, face_tag, param, coords);
    gmsh::model::getDerivative(2, face_tag, param, d1);
    gmsh::model::getSecondDerivative(2, face_tag, param, d2);
    if (coords.size() < 3 || d1.size() < 6 || d2.size() < 9) {
      return false;
    }
    eval.x = {coords[0], coords[1], coords[2]};
    eval.xu = {d1[0], d1[1], d1[2]};
    eval.xv = {d1[3], d1[4], d1[5]};
    eval.xuu = {d2[0], d2[1], d2[2]};
    eval.xvv = {d2[3], d2[4], d2[5]};
    eval.xuv = {d2[6], d2[7], d2[8]};
    return is_finite(eval.x) && is_finite(eval.xu) && is_finite(eval.xv) &&
           is_finite(eval.xuu) && is_finite(eval.xvv) &&
           is_finite(eval.xuv);
  } catch (...) {
    return false;
  }
}

bool surface_uv_for_edge_or_closest(int face_tag, int curve_tag,
                                    double curve_param, const Vec3 &x0,
                                    Vec2 &uv, Vec3 &projected) {
  if (curve_tag > 0 && std::isfinite(curve_param)) {
    try {
      std::vector<double> surface_param;
      gmsh::model::reparametrizeOnSurface(1, curve_tag, {curve_param},
                                          face_tag, surface_param, 0);
      if (surface_param.size() >= 2 && std::isfinite(surface_param[0]) &&
          std::isfinite(surface_param[1])) {
        uv = {surface_param[0], surface_param[1]};
        if (xyz_from_surface_uv(face_tag, uv, projected)) {
          return true;
        }
      }
    } catch (...) {
    }
  }
  return closest_point_on_surface(face_tag, x0, projected, uv);
}

bool surface_normal_cross_norm(int face1, int face2, const Vec3 &x0,
                               int curve_tag, double curve_param,
                               double &cross_norm) {
  Vec3 q1;
  Vec3 q2;
  Vec2 uv1;
  Vec2 uv2;
  if (!surface_uv_for_edge_or_closest(face1, curve_tag, curve_param, x0, uv1,
                                      q1) ||
      !surface_uv_for_edge_or_closest(face2, curve_tag, curve_param, x0, uv2,
                                      q2)) {
    return false;
  }

  SurfaceEval e1;
  SurfaceEval e2;
  if (!surface_eval_d2(face1, uv1, e1) || !surface_eval_d2(face2, uv2, e2)) {
    return false;
  }

  Vec3 n1 = cross(e1.xu, e1.xv);
  Vec3 n2 = cross(e2.xu, e2.xv);
  const double n1_norm = norm(n1);
  const double n2_norm = norm(n2);
  if (n1_norm < 1.0e-14 || n2_norm < 1.0e-14) {
    return false;
  }
  n1 = scale(n1, 1.0 / n1_norm);
  n2 = scale(n2, 1.0 / n2_norm);
  cross_norm = norm(cross(n1, n2));
  return std::isfinite(cross_norm);
}

bool best_face_pair_from_candidates(const std::vector<int> &candidates,
                                    const Vec3 &seed, int &face1, int &face2,
                                    double &normal_cross_norm, int curve_tag,
                                    double curve_param) {
  std::vector<int> faces;
  for (int tag : candidates) {
    append_unique(faces, tag);
  }
  if (faces.size() < 2) {
    return false;
  }

  bool found = false;
  double best_cross = -1.0;
  int best1 = 0;
  int best2 = 0;
  for (std::size_t i = 0; i < faces.size(); ++i) {
    for (std::size_t j = i + 1; j < faces.size(); ++j) {
      double cn = 0.0;
      if (!surface_normal_cross_norm(faces[i], faces[j], seed, curve_tag,
                                     curve_param, cn)) {
        continue;
      }
      if (cn > best_cross) {
        best_cross = cn;
        best1 = faces[i];
        best2 = faces[j];
        found = true;
      }
    }
  }

  if (!found) {
    return false;
  }
  face1 = best1;
  face2 = best2;
  normal_cross_norm = best_cross;
  return true;
}

bool best_adjacent_face_pair(int curve_tag, const Vec3 &seed, int &face1,
                             int &face2, double &normal_cross_norm,
                             const CurveFaceMap *curve_faces = nullptr,
                             double curve_param =
                                 std::numeric_limits<double>::quiet_NaN()) {
  if (curve_faces != nullptr) {
    const auto it = curve_faces->find(curve_tag);
    if (it != curve_faces->end() &&
        best_face_pair_from_candidates(it->second, seed, face1, face2,
                                       normal_cross_norm, curve_tag,
                                       curve_param)) {
      return true;
    }
  }

  std::vector<int> upward;
  std::vector<int> downward;
  try {
    gmsh::model::getAdjacencies(1, curve_tag, upward, downward);
  } catch (...) {
    return false;
  }
  return best_face_pair_from_candidates(upward, seed, face1, face2,
                                        normal_cross_norm, curve_tag,
                                        curve_param);
}

bool solve4x4(double a[4][4], double b[4], double x[4]) {
  double aug[4][5]{};
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      aug[i][j] = a[i][j];
    }
    aug[i][4] = b[i];
  }

  for (int col = 0; col < 4; ++col) {
    int pivot = col;
    double pivot_abs = std::abs(aug[col][col]);
    for (int row = col + 1; row < 4; ++row) {
      const double val = std::abs(aug[row][col]);
      if (val > pivot_abs) {
        pivot = row;
        pivot_abs = val;
      }
    }
    if (pivot_abs < 1.0e-30 || !std::isfinite(pivot_abs)) {
      return false;
    }
    if (pivot != col) {
      for (int j = col; j < 5; ++j) {
        std::swap(aug[col][j], aug[pivot][j]);
      }
    }

    const double diag = aug[col][col];
    for (int j = col; j < 5; ++j) {
      aug[col][j] /= diag;
    }
    for (int row = 0; row < 4; ++row) {
      if (row == col) {
        continue;
      }
      const double factor = aug[row][col];
      for (int j = col; j < 5; ++j) {
        aug[row][j] -= factor * aug[col][j];
      }
    }
  }

  for (int i = 0; i < 4; ++i) {
    x[i] = aug[i][4];
    if (!std::isfinite(x[i])) {
      return false;
    }
  }
  return true;
}

bool edge_newton_4x4_nearest(
    int face1, int face2, const Vec3 &x0, Vec3 &out,
    double tol = 1.0e-14, int maxiter = 25, int curve_tag = -1,
    double curve_param = std::numeric_limits<double>::quiet_NaN(),
    double max_displacement = std::numeric_limits<double>::infinity()) {
  Vec3 p1_proj;
  Vec3 p2_proj;
  Vec2 uv1;
  Vec2 uv2;
  if (!surface_uv_for_edge_or_closest(face1, curve_tag, curve_param, x0, uv1,
                                      p1_proj) ||
      !surface_uv_for_edge_or_closest(face2, curve_tag, curve_param, x0, uv2,
                                      p2_proj)) {
    return false;
  }

  const double match_dist = norm(subtract(p1_proj, p2_proj));
  if (curve_tag > 0 && std::isfinite(curve_param)) {
    return false;
  }
  if (match_dist < std::max(tol * 100.0, 1.0e-5)) {
    out = scale(add(p1_proj, p2_proj), 0.5);
    return true;
  }

  double y[4]{uv1[0], uv1[1], uv2[0], uv2[1]};
  SurfaceEval e1;
  SurfaceEval e2;

  for (int iter = 0; iter < maxiter; ++iter) {
    if (!surface_eval_d2(face1, {y[0], y[1]}, e1) ||
        !surface_eval_d2(face2, {y[2], y[3]}, e2)) {
      return false;
    }

    if (norm(subtract(e1.x, e2.x)) < tol) {
      out = scale(add(e1.x, e2.x), 0.5);
      return true;
    }

    const Vec3 n1 = cross(e1.xu, e1.xv);
    const Vec3 n2 = cross(e2.xu, e2.xv);
    const Vec3 t = cross(n1, n2);

    const Vec3 dn1_du = add(cross(e1.xuu, e1.xv), cross(e1.xu, e1.xuv));
    const Vec3 dn1_dv = add(cross(e1.xuv, e1.xv), cross(e1.xu, e1.xvv));
    const Vec3 dn2_du = add(cross(e2.xuu, e2.xv), cross(e2.xu, e2.xuv));
    const Vec3 dn2_dv = add(cross(e2.xuv, e2.xv), cross(e2.xu, e2.xvv));

    const Vec3 dt_du1 = cross(dn1_du, n2);
    const Vec3 dt_dv1 = cross(dn1_dv, n2);
    const Vec3 dt_du2 = cross(n1, dn2_du);
    const Vec3 dt_dv2 = cross(n1, dn2_dv);

    const Vec3 d = subtract(e1.x, x0);
    const Vec3 xdiff = subtract(e1.x, e2.x);

    double r[4]{xdiff[0], xdiff[1], xdiff[2], dot(d, t)};
    double jac[4][4]{};
    for (int k = 0; k < 3; ++k) {
      jac[k][0] = e1.xu[k];
      jac[k][1] = e1.xv[k];
      jac[k][2] = -e2.xu[k];
      jac[k][3] = -e2.xv[k];
    }
    jac[3][0] = dot(e1.xu, t) + dot(d, dt_du1);
    jac[3][1] = dot(e1.xv, t) + dot(d, dt_dv1);
    jac[3][2] = dot(d, dt_du2);
    jac[3][3] = dot(d, dt_dv2);

    double rhs[4]{-r[0], -r[1], -r[2], -r[3]};
    double dy[4]{};
    if (!solve4x4(jac, rhs, dy)) {
      return false;
    }

    double step_norm = 0.0;
    for (double val : dy) {
      step_norm += val * val;
    }
    step_norm = std::sqrt(step_norm);
    if (step_norm > 2.0) {
      for (double &val : dy) {
        val *= 2.0 / step_norm;
      }
      step_norm = 2.0;
    }

    for (int k = 0; k < 4; ++k) {
      y[k] += dy[k];
      if (!std::isfinite(y[k])) {
        return false;
      }
    }

    double residual_norm = 0.0;
    for (double val : r) {
      residual_norm += val * val;
    }
    residual_norm = std::sqrt(residual_norm);
    if (residual_norm < tol && step_norm < tol) {
      break;
    }
  }

  if (!surface_eval_d2(face1, {y[0], y[1]}, e1) ||
      !surface_eval_d2(face2, {y[2], y[3]}, e2)) {
    return false;
  }
  out = scale(add(e1.x, e2.x), 0.5);
  if (!is_finite(out)) {
    return false;
  }
  if (norm(subtract(e1.x, e2.x)) > 1.0e-8) {
    return false;
  }
  return norm(subtract(out, x0)) <= max_displacement;
}

double get_occ_area() {
  gmsh::vectorpair faces;
  gmsh::model::getEntities(faces, 2);
  double area = 0.0;
  for (const auto &face : faces) {
    try {
      double face_area = 0.0;
      gmsh::model::occ::getMass(2, face.second, face_area);
      if (std::isfinite(face_area)) {
        area += std::abs(face_area);
      }
    } catch (...) {
    }
  }
  return area;
}

double sampled_bbox_diag_from_current_model() {
  double xmin = std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double zmin = std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  double zmax = -std::numeric_limits<double>::infinity();
  bool have_point = false;

  const auto update = [&](const Vec3 &x) {
    if (!is_finite(x)) {
      return;
    }
    xmin = std::min(xmin, x[0]);
    ymin = std::min(ymin, x[1]);
    zmin = std::min(zmin, x[2]);
    xmax = std::max(xmax, x[0]);
    ymax = std::max(ymax, x[1]);
    zmax = std::max(zmax, x[2]);
    have_point = true;
  };

  gmsh::vectorpair curves;
  gmsh::model::getEntities(curves, 1);
  for (const auto &curve : curves) {
    try {
      std::vector<double> pmin;
      std::vector<double> pmax;
      gmsh::model::getParametrizationBounds(1, curve.second, pmin, pmax);
      if (pmin.empty() || pmax.empty() || !std::isfinite(pmin[0]) ||
          !std::isfinite(pmax[0])) {
        continue;
      }
      constexpr int n_curve = 64;
      for (int i = 0; i <= n_curve; ++i) {
        const double a = static_cast<double>(i) / n_curve;
        const double t = (1.0 - a) * pmin[0] + a * pmax[0];
        std::vector<double> coords;
        gmsh::model::getValue(1, curve.second, {t}, coords);
        if (coords.size() >= 3) {
          update({coords[0], coords[1], coords[2]});
        }
      }
    } catch (...) {
    }
  }

  gmsh::vectorpair faces;
  gmsh::model::getEntities(faces, 2);
  for (const auto &face : faces) {
    try {
      std::vector<double> pmin;
      std::vector<double> pmax;
      gmsh::model::getParametrizationBounds(2, face.second, pmin, pmax);
      if (pmin.size() < 2 || pmax.size() < 2 || !std::isfinite(pmin[0]) ||
          !std::isfinite(pmin[1]) || !std::isfinite(pmax[0]) ||
          !std::isfinite(pmax[1])) {
        continue;
      }
      constexpr int n_face = 16;
      for (int i = 0; i <= n_face; ++i) {
        const double ai = static_cast<double>(i) / n_face;
        const double u = (1.0 - ai) * pmin[0] + ai * pmax[0];
        for (int j = 0; j <= n_face; ++j) {
          const double aj = static_cast<double>(j) / n_face;
          const double v = (1.0 - aj) * pmin[1] + aj * pmax[1];
          std::vector<double> coords;
          gmsh::model::getValue(2, face.second, {u, v}, coords);
          if (coords.size() >= 3) {
            update({coords[0], coords[1], coords[2]});
          }
        }
      }
    } catch (...) {
    }
  }

  if (!have_point) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  constexpr double cadquery_bbox_padding = 9.99843e-8;
  const double dx = (xmax - xmin) + 2.0 * cadquery_bbox_padding;
  const double dy = (ymax - ymin) + 2.0 * cadquery_bbox_padding;
  const double dz = (zmax - zmin) + 2.0 * cadquery_bbox_padding;
  const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
  return std::isfinite(diag) && diag > 0.0
             ? diag
             : std::numeric_limits<double>::quiet_NaN();
}

void append_xyz(std::vector<double> &data, const Vec3 &xyz) {
  data.push_back(xyz[0]);
  data.push_back(xyz[1]);
  data.push_back(xyz[2]);
}

CurveFaceMap build_curve_face_map_from_mesh(
    const std::unordered_map<int, int> &current_face_to_original_face) {
  std::unordered_map<EdgeKey, std::vector<int>, EdgeKeyHash> edge_to_faces;

  gmsh::vectorpair faces;
  gmsh::model::getEntities(faces, 2);
  for (const auto &face : faces) {
    const int face_tag = face.second;
    const auto original_it = current_face_to_original_face.find(face_tag);
    if (original_it == current_face_to_original_face.end()) {
      continue;
    }
    const int original_face_tag = original_it->second;
    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;
    gmsh::model::mesh::getElements(element_types, element_tags,
                                   element_node_tags, 2, face_tag);

    for (std::size_t block = 0; block < element_types.size(); ++block) {
      std::string element_name;
      int dim = 0;
      int element_order = 0;
      int num_nodes = 0;
      int num_primary_nodes = 0;
      std::vector<double> local_coords;
      gmsh::model::mesh::getElementProperties(
          element_types[block], element_name, dim, element_order, num_nodes,
          local_coords, num_primary_nodes);
      if (dim != 2 || num_primary_nodes < 3 ||
          element_name.find("Triangle") == std::string::npos) {
        continue;
      }

      const auto &tags = element_node_tags[block];
      const std::size_t num_elements =
          tags.size() / static_cast<std::size_t>(num_nodes);
      for (std::size_t elem = 0; elem < num_elements; ++elem) {
        const auto tag0 = tags[elem * static_cast<std::size_t>(num_nodes)];
        const auto tag1 =
            tags[elem * static_cast<std::size_t>(num_nodes) + 1];
        const auto tag2 =
            tags[elem * static_cast<std::size_t>(num_nodes) + 2];
        append_unique(edge_to_faces[unoriented_edge(tag0, tag1)],
                      original_face_tag);
        append_unique(edge_to_faces[unoriented_edge(tag1, tag2)],
                      original_face_tag);
        append_unique(edge_to_faces[unoriented_edge(tag2, tag0)],
                      original_face_tag);
      }
    }
  }

  CurveFaceMap curve_faces;
  gmsh::vectorpair curves;
  gmsh::model::getEntities(curves, 1);
  for (const auto &curve : curves) {
    const int curve_tag = curve.second;
    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;
    gmsh::model::mesh::getElements(element_types, element_tags,
                                   element_node_tags, 1, curve_tag);

    for (std::size_t block = 0; block < element_types.size(); ++block) {
      std::string element_name;
      int dim = 0;
      int element_order = 0;
      int num_nodes = 0;
      int num_primary_nodes = 0;
      std::vector<double> local_coords;
      gmsh::model::mesh::getElementProperties(
          element_types[block], element_name, dim, element_order, num_nodes,
          local_coords, num_primary_nodes);
      if (dim != 1 || num_primary_nodes < 2 ||
          element_name.find("Line") == std::string::npos) {
        continue;
      }

      const auto &tags = element_node_tags[block];
      const std::size_t num_elements =
          tags.size() / static_cast<std::size_t>(num_nodes);
      for (std::size_t elem = 0; elem < num_elements; ++elem) {
        const auto tag0 = tags[elem * static_cast<std::size_t>(num_nodes)];
        const auto tag1 =
            tags[elem * static_cast<std::size_t>(num_nodes) + 1];
        const auto faces_it = edge_to_faces.find(unoriented_edge(tag0, tag1));
        if (faces_it == edge_to_faces.end()) {
          continue;
        }
        for (int face_tag : faces_it->second) {
          append_unique(curve_faces[curve_tag], face_tag);
        }
      }
    }
  }

  return curve_faces;
}

struct GllEdgePanel {
  std::size_t tag0 = 0;
  std::size_t tag1 = 0;
  double tmid = std::numeric_limits<double>::quiet_NaN();
  std::vector<Vec3> samples;
};

double curve_node_param(std::size_t tag) {
  try {
    std::vector<double> coord;
    std::vector<double> param;
    int dim = -1;
    int ent = -1;
    gmsh::model::mesh::getNode(tag, coord, param, dim, ent);
    if (!param.empty() && std::isfinite(param[0])) {
      return param[0];
    }
  } catch (...) {
  }
  return std::numeric_limits<double>::quiet_NaN();
}

void append_panel_data(
    MeshResult &result, std::size_t tag0, std::size_t tag1,
    const std::vector<double> &edge_nodes,
    const std::vector<double> &interp_weights,
    const std::vector<double> &quad_weights, const std::vector<Vec3> &samples) {
  if (samples.size() != edge_nodes.size() || samples.size() != quad_weights.size()) {
    return;
  }

  const auto derivatives =
      edge_derivatives_at_nodes(edge_nodes, interp_weights, samples);

  result.panel_node_tags.push_back(static_cast<int>(tag0));
  result.panel_node_tags.push_back(static_cast<int>(tag1));

  for (std::size_t i = 0; i < samples.size(); ++i) {
    Vec3 d = derivatives[i];
    double speed = norm(d);
    if (!std::isfinite(speed) || speed <= 0.0) {
      d = finite_difference_edge_tangent(edge_nodes, samples, i);
      speed = norm(d);
    }

    Vec3 tangent{1.0, 0.0, 0.0};
    double line_weight = 0.0;
    if (std::isfinite(speed) && speed > 0.0) {
      tangent = scale(d, 1.0 / speed);
      line_weight = quad_weights[i] * speed;
    }

    append_xyz(result.panel_xyz, samples[i]);
    append_xyz(result.panel_tangents, tangent);
    result.panel_weights.push_back(line_weight);
  }

  ++result.npanels_total;
}

void record_edge_panels(
    const std::vector<double> &edge_nodes,
    const std::vector<double> &interp_weights,
    const std::vector<double> &quad_weights, const EdgeSampleMap &edge_samples,
    const CurveFaceMap &curve_faces, const FeatureCurveMap &feature_curves,
    MeshResult &result) {
  result.edge_order = static_cast<int>(edge_nodes.size()) - 1;
  result.edge_nodes_per_panel = static_cast<int>(edge_nodes.size());
  result.edge_gll_nodes = edge_nodes;

  gmsh::vectorpair curves;
  gmsh::model::getEntities(curves, 1);
  for (const auto &curve : curves) {
    const int curve_tag = curve.second;
    if (feature_curves.find(curve_tag) != feature_curves.end()) {
      continue;
    }

    std::vector<int> faces;
    const auto faces_it = curve_faces.find(curve_tag);
    if (faces_it != curve_faces.end()) {
      faces = faces_it->second;
    }

    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;
    gmsh::model::mesh::getElements(element_types, element_tags,
                                   element_node_tags, 1, curve_tag);

    std::vector<GllEdgePanel> panels;
    for (std::size_t block = 0; block < element_types.size(); ++block) {
      std::string element_name;
      int dim = 0;
      int element_order = 0;
      int num_nodes = 0;
      int num_primary_nodes = 0;
      std::vector<double> local_coords;
      gmsh::model::mesh::getElementProperties(
          element_types[block], element_name, dim, element_order, num_nodes,
          local_coords, num_primary_nodes);
      if (dim != 1 || num_primary_nodes < 2 ||
          element_name.find("Line") == std::string::npos) {
        continue;
      }

      const auto &tags = element_node_tags[block];
      const std::size_t num_elements =
          tags.size() / static_cast<std::size_t>(num_nodes);
      for (std::size_t elem = 0; elem < num_elements; ++elem) {
        const std::size_t tag0 =
            tags[elem * static_cast<std::size_t>(num_nodes)];
        const std::size_t tag1 =
            tags[elem * static_cast<std::size_t>(num_nodes) + 1];
        std::size_t out_tag0 = tag0;
        std::size_t out_tag1 = tag1;
        double t0 = curve_node_param(tag0);
        double t1 = curve_node_param(tag1);
        if (std::isfinite(t0) && std::isfinite(t1) && t1 < t0) {
          std::swap(out_tag0, out_tag1);
          std::swap(t0, t1);
        }

        const auto edge_it = edge_samples.find(EdgeKey{out_tag0, out_tag1});
        if (edge_it == edge_samples.end()) {
          continue;
        }
        GllEdgePanel panel;
        panel.tag0 = out_tag0;
        panel.tag1 = out_tag1;
        if (std::isfinite(t0) && std::isfinite(t1)) {
          panel.tmid = 0.5 * (t0 + t1);
        }
        panel.samples = edge_it->second;
        panels.push_back(std::move(panel));
      }
    }

    if (panels.empty()) {
      continue;
    }

    std::stable_sort(panels.begin(), panels.end(),
                     [](const GllEdgePanel &a, const GllEdgePanel &b) {
                       if (std::isfinite(a.tmid) && std::isfinite(b.tmid)) {
                         return a.tmid < b.tmid;
                       }
                       if (std::isfinite(a.tmid) != std::isfinite(b.tmid)) {
                         return std::isfinite(a.tmid);
                       }
                       return a.tag0 < b.tag0;
                     });

    std::vector<int> adjacent{-1, -1};
    for (std::size_t i = 0; i < faces.size() && i < 2; ++i) {
      adjacent[i] = faces[i];
    }

    ++result.ncad_edges;
    result.edge_ids.push_back(result.ncad_edges);
    result.edge_curve_tags.push_back(curve_tag);
    result.edge_adjacent_faces.push_back(adjacent[0]);
    result.edge_adjacent_faces.push_back(adjacent[1]);
    result.edge_panel_start.push_back(result.npanels_total + 1);
    result.edge_panel_count.push_back(static_cast<int>(panels.size()));

    for (const GllEdgePanel &panel : panels) {
      append_panel_data(result, panel.tag0, panel.tag1, edge_nodes,
                        interp_weights, quad_weights, panel.samples);
    }
  }
}

FeatureNodeConstraintMap
collect_feature_node_constraints(const FeatureCurveMap &feature_curves) {
  FeatureNodeConstraintMap constraints;
  for (const auto &entry : feature_curves) {
    const int curve_tag = entry.first;
    const FeatureCurveInfo &info = entry.second;

    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;
    gmsh::model::mesh::getElements(element_types, element_tags,
                                   element_node_tags, 1, curve_tag);

    for (std::size_t block = 0; block < element_types.size(); ++block) {
      std::string element_name;
      int dim = 0;
      int element_order = 0;
      int num_nodes = 0;
      int num_primary_nodes = 0;
      std::vector<double> local_coords;
      gmsh::model::mesh::getElementProperties(
          element_types[block], element_name, dim, element_order, num_nodes,
          local_coords, num_primary_nodes);
      if (dim != 1 || num_primary_nodes < 2 ||
          element_name.find("Line") == std::string::npos) {
        continue;
      }

      const auto &tags = element_node_tags[block];
      const std::size_t num_elements =
          tags.size() / static_cast<std::size_t>(num_nodes);
      for (std::size_t elem = 0; elem < num_elements; ++elem) {
        for (int local = 0; local < num_primary_nodes; ++local) {
          const std::size_t tag =
              tags[elem * static_cast<std::size_t>(num_nodes) +
                   static_cast<std::size_t>(local)];
          auto &node_constraints = constraints[tag];
          const bool already_present = std::any_of(
              node_constraints.begin(), node_constraints.end(),
              [&](const FeatureNodeConstraint &constraint) {
                return constraint.original_face_tag == info.original_face_tag &&
                       constraint.axis == info.axis &&
                       std::abs(constraint.param - info.param) <=
                           1.0e-13 *
                               std::max({1.0, std::abs(constraint.param),
                                         std::abs(info.param)});
              });
          if (!already_present) {
            node_constraints.push_back(
                {info.original_face_tag, info.axis, info.param});
          }
        }
      }
    }
  }
  return constraints;
}

bool feature_node_should_use_constraint(std::size_t node_tag, int dim,
                                        int entity_tag,
                                        const FeatureCurveMap &feature_curves,
                                        const FeatureNodeConstraintMap
                                            &constraints) {
  if (constraints.find(node_tag) == constraints.end()) {
    return false;
  }
  if (dim == 2) {
    return true;
  }
  if (dim == 1) {
    return feature_curves.find(entity_tag) != feature_curves.end();
  }
  if (dim == 0) {
    std::vector<int> upward;
    std::vector<int> downward;
    try {
      gmsh::model::getAdjacencies(0, entity_tag, upward, downward);
    } catch (...) {
      return false;
    }
    bool has_feature_curve = false;
    bool has_true_curve = false;
    for (int curve_tag : upward) {
      if (feature_curves.find(curve_tag) != feature_curves.end()) {
        has_feature_curve = true;
      } else {
        has_true_curve = true;
      }
    }
    return has_feature_curve && !has_true_curve;
  }
  return false;
}

std::unordered_map<std::size_t, NodeRecord>
get_node_map(MeshResult &result, const CurveFaceMap &curve_faces,
             const FeatureCurveMap &feature_curves,
             occ::Context &cad, bool use_occ) {
  std::vector<std::size_t> node_tags;
  std::vector<double> coords;
  std::vector<double> params;
  gmsh::model::mesh::getNodes(node_tags, coords, params, -1, -1, false, false);

  std::unordered_map<std::size_t, Vec3> raw_nodes;
  raw_nodes.reserve(node_tags.size());
  for (std::size_t i = 0; i < node_tags.size(); ++i) {
    raw_nodes[node_tags[i]] =
        Vec3{coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]};
  }
  const FeatureNodeConstraintMap feature_node_constraints =
      collect_feature_node_constraints(feature_curves);
  const auto feature_node_helper = load_occ_feature_nodes(
      cad, use_occ, raw_nodes, feature_node_constraints);

  std::unordered_map<std::size_t, NodeRecord> nodes;
  nodes.reserve(node_tags.size());
  result.nscaffold_nodes = static_cast<int>(node_tags.size());
  result.scaffold_nodes.reserve(3 * node_tags.size());
  result.scaffold_node_tags.reserve(node_tags.size());

  for (std::size_t i = 0; i < node_tags.size(); ++i) {
    Vec3 xyz{coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]};
    Vec3 corrected = xyz;

    try {
      std::vector<double> node_coord;
      std::vector<double> node_param;
      int dim = -1;
      int entity_tag = -1;
      gmsh::model::mesh::getNode(node_tags[i], node_coord, node_param, dim,
                                 entity_tag);
      if (node_coord.size() >= 3 && std::isfinite(node_coord[0]) &&
          std::isfinite(node_coord[1]) && std::isfinite(node_coord[2])) {
        xyz = {node_coord[0], node_coord[1], node_coord[2]};
        corrected = xyz;
      }

      const bool use_feature_constraint = feature_node_helper.ok &&
          feature_node_should_use_constraint(node_tags[i], dim, entity_tag,
                                             feature_curves,
                                             feature_node_constraints);
      const auto feature_node_it =
          use_feature_constraint ? feature_node_helper.nodes.find(node_tags[i])
                                 : feature_node_helper.nodes.end();
      if (feature_node_it != feature_node_helper.nodes.end()) {
        corrected = feature_node_it->second;
      } else if (dim == 2) {
        Vec2 uv{};
        Vec3 projected{};
        if (closest_point_on_surface(entity_tag, xyz, projected, uv)) {
          corrected = projected;
        }
      } else if (dim == 1 &&
                 feature_curves.find(entity_tag) == feature_curves.end()) {
        Vec3 seed_on_curve = xyz;
        double curve_param =
            !node_param.empty() && std::isfinite(node_param[0])
                ? node_param[0]
                : std::numeric_limits<double>::quiet_NaN();
        if (closest_point_on_curve(entity_tag, xyz, seed_on_curve,
                                   &curve_param)) {
          corrected = seed_on_curve;
          ++result.step1b_edge_curve_projection_success;
        }

        int face1 = 0;
        int face2 = 0;
        double normal_cross = 0.0;
        if (best_adjacent_face_pair(entity_tag, seed_on_curve, face1, face2,
                                    normal_cross, &curve_faces,
                                    curve_param) &&
            normal_cross > 1.0e-8) {
          Vec3 newton_point{};
          const double max_displacement =
              std::max(1.0e-8, result.mesh_size);
          if (edge_newton_4x4_nearest(face1, face2, seed_on_curve,
                                      newton_point, 1.0e-14, 25,
                                      entity_tag, curve_param,
                                      max_displacement)) {
            corrected = newton_point;
            ++result.step1b_edge_newton_success;
          }
        }
      }
    } catch (...) {
      corrected = xyz;
    }

    nodes.emplace(node_tags[i],
                  NodeRecord{corrected, static_cast<int>(i) + 1});
    result.scaffold_node_tags.push_back(static_cast<int>(node_tags[i]));
    append_xyz(result.scaffold_nodes, corrected);
  }
  return nodes;
}

EdgeSampleMap build_edge_samples(
    const std::unordered_map<std::size_t, NodeRecord> &node_map,
    const std::vector<double> &edge_nodes, const CurveFaceMap &curve_faces,
    const FeatureCurveMap &feature_curves, MeshResult &result) {
  EdgeSampleMap edge_samples;
  gmsh::vectorpair curves;
  gmsh::model::getEntities(curves, 1);

  for (const auto &curve : curves) {
    const int curve_tag = curve.second;
    if (feature_curves.find(curve_tag) != feature_curves.end()) {
      continue;
    }
    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;
    gmsh::model::mesh::getElements(element_types, element_tags,
                                   element_node_tags, 1, curve_tag);

    for (std::size_t block = 0; block < element_types.size(); ++block) {
      std::string element_name;
      int dim = 0;
      int element_order = 0;
      int num_nodes = 0;
      int num_primary_nodes = 0;
      std::vector<double> local_coords;
      gmsh::model::mesh::getElementProperties(
          element_types[block], element_name, dim, element_order, num_nodes,
          local_coords, num_primary_nodes);

      if (dim != 1 || num_primary_nodes < 2 ||
          element_name.find("Line") == std::string::npos) {
        continue;
      }

      const auto &tags = element_node_tags[block];
      const std::size_t num_elements =
          tags.size() / static_cast<std::size_t>(num_nodes);
      for (std::size_t elem = 0; elem < num_elements; ++elem) {
        const auto tag0 = tags[elem * static_cast<std::size_t>(num_nodes)];
        const auto tag1 =
            tags[elem * static_cast<std::size_t>(num_nodes) + 1];
        const auto p0_it = node_map.find(tag0);
        const auto p1_it = node_map.find(tag1);
        if (p0_it == node_map.end() || p1_it == node_map.end()) {
          continue;
        }

        const Vec3 p0 = p0_it->second.xyz;
        const Vec3 p1 = p1_it->second.xyz;
        double t0 = std::numeric_limits<double>::quiet_NaN();
        double t1 = std::numeric_limits<double>::quiet_NaN();
        try {
          std::vector<double> coord0;
          std::vector<double> coord1;
          std::vector<double> param0;
          std::vector<double> param1;
          int dim0 = -1;
          int dim1 = -1;
          int ent0 = -1;
          int ent1 = -1;
          gmsh::model::mesh::getNode(tag0, coord0, param0, dim0, ent0);
          gmsh::model::mesh::getNode(tag1, coord1, param1, dim1, ent1);
          if (!param0.empty() && !param1.empty() && std::isfinite(param0[0]) &&
              std::isfinite(param1[0])) {
            t0 = param0[0];
            t1 = param1[0];
          }
        } catch (...) {
        }

        std::vector<Vec3> samples;
        samples.reserve(edge_nodes.size());
        for (const double a : edge_nodes) {
          if (std::abs(a) < 1.0e-15) {
            samples.push_back(p0);
            continue;
          }
          if (std::abs(a - 1.0) < 1.0e-15) {
            samples.push_back(p1);
            continue;
          }

          const Vec3 linear = blend_segment(p0, p1, a);
          Vec3 seed_on_curve = linear;
          double curve_param =
              std::isfinite(t0) && std::isfinite(t1)
                  ? (1.0 - a) * t0 + a * t1
                  : std::numeric_limits<double>::quiet_NaN();
          bool used_curve_projection = false;
          if (closest_point_on_curve(curve_tag, linear, seed_on_curve,
                                     &curve_param)) {
            used_curve_projection = true;
            ++result.edge_curve_projection_success;
          } else {
            ++result.edge_fallback_points;
          }

          Vec3 final_point = seed_on_curve;
          int face1 = 0;
          int face2 = 0;
          double normal_cross = 0.0;
          if (best_adjacent_face_pair(curve_tag, seed_on_curve, face1, face2,
                                      normal_cross, &curve_faces,
                                      curve_param) &&
              normal_cross > 1.0e-8) {
            Vec3 newton_point{};
            const double max_displacement =
                std::max(1.0e-8, 0.75 * norm(subtract(p1, p0)));
            if (edge_newton_4x4_nearest(face1, face2, seed_on_curve,
                                        newton_point, 1.0e-14, 25,
                                        curve_tag, curve_param,
                                        max_displacement)) {
              final_point = newton_point;
              ++result.edge_newton_success;
            }
          }

          if (!used_curve_projection && !is_finite(final_point)) {
            final_point = linear;
          }
          samples.push_back(final_point);
        }

        edge_samples[EdgeKey{tag0, tag1}] = samples;
        std::vector<Vec3> reversed = samples;
        std::reverse(reversed.begin(), reversed.end());
        edge_samples[EdgeKey{tag1, tag0}] = reversed;
      }
    }
  }

  return edge_samples;
}

int add_helper_curve_to_gmsh(const std::vector<Vec3> &raw_points) {
  std::vector<Vec3> points;
  points.reserve(raw_points.size());
  constexpr double point_tol = 1.0e-11;
  for (const Vec3 &p : raw_points) {
    if (points.empty() || norm(subtract(p, points.back())) > point_tol) {
      points.push_back(p);
    }
  }
  if (points.size() < 2) {
    return -1;
  }
  if (norm(subtract(points.front(), points.back())) <= point_tol &&
      points.size() >= 2) {
    points.pop_back();
    if (points.size() < 3) {
      return -1;
    }
  }

  double length = 0.0;
  for (std::size_t i = 1; i < points.size(); ++i) {
    length += norm(subtract(points[i], points[i - 1]));
  }
  if (length <= 1.0e-12) {
    return -1;
  }

  std::vector<int> point_tags;
  point_tags.reserve(points.size());
  for (const Vec3 &p : points) {
    point_tags.push_back(gmsh::model::occ::addPoint(p[0], p[1], p[2], 0.0));
  }
  if (point_tags.size() == 2) {
    return gmsh::model::occ::addLine(point_tags[0], point_tags[1]);
  }
  return gmsh::model::occ::addSpline(point_tags);
}

FragmentationInfo apply_continuity_fragmentation(
    occ::Context &cad, bool use_occ,
    double mesh_size, const std::vector<int> &original_face_tags) {
  FragmentationInfo info;
  for (int face_tag : original_face_tags) {
    info.current_face_to_original_face[face_tag] = face_tag;
  }

  const double helper_spacing = std::max(1.0e-10, 0.75 * mesh_size);
  auto helper = load_occ_helper_curves(
      cad, use_occ, helper_spacing, original_face_tags);
  if (!helper.ok || helper.curves.empty()) {
    if (helper.ok) {
      info.continuity_break_count = helper.break_count;
    }
    return info;
  }

  info.continuity_break_count = helper.break_count;

  std::vector<int> helper_curve_tags;
  std::vector<FeatureCurveInfo> helper_infos;
  helper_curve_tags.reserve(helper.curves.size());
  helper_infos.reserve(helper.curves.size());
  for (const ContinuityHelperCurve &curve : helper.curves) {
    const int curve_tag = add_helper_curve_to_gmsh(curve.points);
    if (curve_tag <= 0) {
      continue;
    }
    helper_curve_tags.push_back(curve_tag);
    helper_infos.push_back(
        {curve.original_face_tag, curve.axis, curve.param});
  }
  info.helper_curve_count = static_cast<int>(helper_curve_tags.size());
  if (helper_curve_tags.empty()) {
    return info;
  }

  gmsh::model::occ::synchronize();
  gmsh::vectorpair object_faces;
  object_faces.reserve(original_face_tags.size());
  for (int face_tag : original_face_tags) {
    object_faces.push_back({2, face_tag});
  }
  gmsh::vectorpair tool_curves;
  tool_curves.reserve(helper_curve_tags.size());
  for (int curve_tag : helper_curve_tags) {
    tool_curves.push_back({1, curve_tag});
  }

  gmsh::vectorpair out_dim_tags;
  std::vector<gmsh::vectorpair> out_dim_tags_map;
  gmsh::model::occ::fragment(object_faces, tool_curves, out_dim_tags,
                             out_dim_tags_map, -1, true, true);
  gmsh::model::occ::synchronize();

  info.current_face_to_original_face.clear();
  for (std::size_t i = 0; i < original_face_tags.size(); ++i) {
    if (i >= out_dim_tags_map.size()) {
      continue;
    }
    std::vector<int> fragment_faces;
    for (const auto &dim_tag : out_dim_tags_map[i]) {
      if (dim_tag.first == 2) {
        fragment_faces.push_back(dim_tag.second);
        const auto inserted = info.current_face_to_original_face.emplace(
            dim_tag.second, original_face_tags[i]);
        if (!inserted.second && inserted.first->second != original_face_tags[i])
          throw std::runtime_error("[native_identity] Fragment has distinct original parents: face " + std::to_string(dim_tag.second));
      }
    }
    if (fragment_faces.size() != 1 ||
        fragment_faces.front() != original_face_tags[i]) {
      ++info.fragmented_face_count;
    }
  }

  const std::size_t helper_offset = object_faces.size();
  for (std::size_t i = 0; i < helper_curve_tags.size(); ++i) {
    const std::size_t map_idx = helper_offset + i;
    if (map_idx >= out_dim_tags_map.size()) {
      continue;
    }
    for (const auto &dim_tag : out_dim_tags_map[map_idx]) {
      if (dim_tag.first == 1) {
        info.feature_curves[dim_tag.second] = helper_infos[i];
      }
    }
  }

  for (const auto &entry : info.feature_curves) {
    const int curve_tag = entry.first;
    const int original_face = entry.second.original_face_tag;
    std::vector<int> target_faces;
    for (const auto &face_entry : info.current_face_to_original_face) {
      if (face_entry.second == original_face) {
        target_faces.push_back(face_entry.first);
      }
    }

    std::vector<int> upward;
    std::vector<int> downward;
    try {
      gmsh::model::getAdjacencies(1, curve_tag, upward, downward);
    } catch (...) {
      upward.clear();
    }
    bool already_on_target = false;
    for (int face_tag : upward) {
      const auto it = info.current_face_to_original_face.find(face_tag);
      if (it != info.current_face_to_original_face.end() &&
          it->second == original_face) {
        already_on_target = true;
        break;
      }
    }
    if (already_on_target) {
      continue;
    }
    for (int face_tag : target_faces) {
      try {
        gmsh::model::mesh::embed(1, {curve_tag}, 2, face_tag);
      } catch (...) {
      }
    }
  }

  return info;
}

double gmsh_algorithm_id(const std::string &name) {
  if (name == "frontal_delaunay") {
    return 6.0;
  }
  if (name == "delaunay") {
    return 5.0;
  }
  if (name == "meshadapt") {
    return 1.0;
  }
  if (name == "automatic") {
    return 2.0;
  }
  throw std::invalid_argument("Unknown gmsh_algorithm: " + name);
}

void validate_mesh_options(const MeshOptions &options) {
  if (options.mesh_profile == "curvature") {
    if (options.curvature_points <= 0) {
      throw std::invalid_argument("curvature_points must be positive");
    }
    if (options.hmin_mode != "fraction" && options.hmin_mode != "cad_edge") {
      throw std::invalid_argument("hmin_mode must be 'fraction' or 'cad_edge'");
    }
    if (!(options.hmin_fraction > 0.0) || !(options.hmax_fraction > 0.0)) {
      throw std::invalid_argument("hmin_fraction and hmax_fraction must be positive");
    }
    if (!(options.hmin_edge_scale > 0.0)) {
      throw std::invalid_argument("hmin_edge_scale must be positive");
    }
    if (!(options.cad_edge_tol_fraction >= 0.0)) {
      throw std::invalid_argument("cad_edge_tol_fraction must be nonnegative");
    }
    if (options.hmin_fraction > options.hmax_fraction) {
      throw std::invalid_argument("hmin_fraction must be <= hmax_fraction");
    }
    (void)gmsh_algorithm_id(options.gmsh_algorithm);
  }
  if (!(options.min_angle_deg > 0.0) || !(options.min_angle_deg < 60.0)) {
    throw std::invalid_argument("min_angle_deg must be between 0 and 60");
  }
  if (!(options.max_edge_ratio >= 1.0)) {
    throw std::invalid_argument("max_edge_ratio must be >= 1");
  }
}

double effective_hmin(double bbox_diag, const MeshOptions &options,
                      const OccInfo &occ) {
  if (options.hmin_mode == "fraction") {
    return options.hmin_fraction * bbox_diag;
  }
  if (!occ.ok ||
      !std::isfinite(occ.min_meaningful_cad_edge_length) ||
      occ.min_meaningful_cad_edge_length <= 0.0 ||
      occ.meaningful_cad_edge_count <= 0) {
    throw std::runtime_error(
        "hmin_mode='cad_edge' requested, but no meaningful CAD edge length "
        "was available from OCC. Try hmin_mode='fraction' or relax "
        "cad_edge_tol_fraction.");
  }
  return options.hmin_edge_scale * occ.min_meaningful_cad_edge_length;
}

void configure_meshing(double mesh_size, double hmin, double hmax,
                       const MeshOptions &options) {
  set_option_if_available("Mesh.CharacteristicLengthMin", mesh_size);
  set_option_if_available("Mesh.CharacteristicLengthMax", mesh_size);
  set_option_if_available("Mesh.ElementOrder", 1.0);
  set_option_if_available("Mesh.SecondOrderLinear", 0.0);
  set_option_if_available("Mesh.RecombineAll", 0.0);
  set_option_if_available("Mesh.MeshSizeFromPoints", 0.0);
  set_option_if_available("Mesh.MeshSizeFromCurvature", 0.0);
  set_option_if_available("Mesh.MeshSizeExtendFromBoundary", 0.0);

  if (options.mesh_profile == "quality") {
    set_option_if_available("Mesh.Algorithm", 6.0);
    set_option_if_available("Mesh.Optimize", 1.0);
  } else if (options.mesh_profile == "adaptive") {
    set_option_if_available("Mesh.Algorithm", 5.0);
  } else if (options.mesh_profile == "robust") {
    set_option_if_available("Mesh.Algorithm", 2.0);
  } else if (options.mesh_profile == "curvature") {
    set_option_if_available("Mesh.CharacteristicLengthMin", hmin);
    set_option_if_available("Mesh.CharacteristicLengthMax", hmax);
    set_option_if_available("Mesh.MeshSizeFromPoints", 1.0);
    set_option_if_available("Mesh.MeshSizeFromCurvature",
                            static_cast<double>(options.curvature_points));
    set_option_if_available("Mesh.MeshSizeExtendFromBoundary",
                            options.mesh_size_extend_from_boundary ? 1.0
                                                                   : 0.0);
    set_option_if_available("Mesh.Algorithm",
                            gmsh_algorithm_id(options.gmsh_algorithm));
    set_option_if_available("Mesh.Optimize", options.optimize ? 1.0 : 0.0);
  }
}

} // namespace

MeshResult mesh_step_file(const std::string &step_file,
                          const MeshOptions &options,
                          const SampleNodes &sample_nodes) {
  if (sample_nodes.nodes_per_patch <= 0 ||
      sample_nodes.uv.size() !=
          2 * static_cast<std::size_t>(sample_nodes.nodes_per_patch)) {
    throw std::invalid_argument("sample_nodes must contain a 2 x N uv array");
  }
  validate_mesh_options(options);

  std::lock_guard<std::mutex> cad_lock(cad_process_mutex);
  // Destroy Gmsh bindings before releasing the retained shape.
  std::unique_ptr<occ::Context> retained;
  GmshSession session(options.verbose);
  verify_native_runtime();
  const bool use_occ = use_occ_computations(options.occ_info_exe);
  retained = std::make_unique<occ::Context>(step_file);
  auto &cad = *retained;
  MeshResult result;
  result.step_file = step_file;
  result.order = options.order;
  result.nodes_per_patch = sample_nodes.nodes_per_patch;
  result.mesh_profile = options.mesh_profile;
  result.sample_rule = "caller_provided";

  const OccInfo occ = load_occ_info(cad, use_occ, options.cad_edge_tol_fraction);

  gmsh::model::add("step_mesher_cpp_matlab");
  const auto bound_face_tags = cad.bind_native();
  result.native_identity_faces = static_cast<int>(bound_face_tags.size());
  result.native_identity_occurrences = static_cast<int>(cad.face_occurrences().size());
  result.occ_evaluation = use_occ ? "original_occ_in_process" : "original_gmsh";
  gmsh::model::occ::synchronize();

  gmsh::vectorpair volumes;
  gmsh::model::getEntities(volumes, 3);
  if (!volumes.empty()) {
    gmsh::model::occ::remove(volumes, false);
    gmsh::model::occ::synchronize();
  }

  gmsh::vectorpair faces;
  gmsh::model::getEntities(faces, 2);
  if (faces.empty()) {
    throw std::runtime_error("No CAD faces were imported from STEP file");
  }
  std::vector<int> original_face_tags;
  original_face_tags.reserve(faces.size());
  for (const auto &face : faces) {
    original_face_tags.push_back(face.second);
  }
  if (original_face_tags.size() != bound_face_tags.size())
    throw std::runtime_error("[native_identity] Surface coverage changed after volume removal");

  double xmin = 0.0;
  double ymin = 0.0;
  double zmin = 0.0;
  double xmax = 0.0;
  double ymax = 0.0;
  double zmax = 0.0;
  gmsh::model::getBoundingBox(-1, -1, xmin, ymin, zmin, xmax, ymax, zmax);
  const double sampled_bbox_diag = sampled_bbox_diag_from_current_model();
  const double gmsh_bbox_diag =
      std::sqrt((xmax - xmin) * (xmax - xmin) +
                (ymax - ymin) * (ymax - ymin) +
                (zmax - zmin) * (zmax - zmin));
  result.bbox_diag = occ.ok && std::isfinite(occ.tight_bbox_diag) &&
                             occ.tight_bbox_diag > 0.0
                         ? occ.tight_bbox_diag
                         : (std::isfinite(sampled_bbox_diag) &&
                                    sampled_bbox_diag > 0.0
                                ? sampled_bbox_diag
                                : gmsh_bbox_diag);
  result.mesh_size = std::max(options.mesh_fraction * result.bbox_diag,
                              std::numeric_limits<double>::epsilon());
  result.min_raw_cad_edge_length = occ.min_raw_cad_edge_length;
  result.min_meaningful_cad_edge_length =
      occ.min_meaningful_cad_edge_length;
  result.cad_edge_count = occ.cad_edge_count;
  result.meaningful_cad_edge_count = occ.meaningful_cad_edge_count;
  result.hmin_effective = result.mesh_size;
  result.hmax_effective = result.mesh_size;
  if (options.mesh_profile == "curvature") {
    result.hmin_effective = effective_hmin(result.bbox_diag, options, occ);
    result.hmax_effective = options.hmax_fraction * result.bbox_diag;
    if (!(result.hmin_effective > 0.0) ||
        !(result.hmax_effective > 0.0) ||
        !std::isfinite(result.hmin_effective) ||
        !std::isfinite(result.hmax_effective)) {
      throw std::runtime_error("computed non-finite curvature mesh size bound");
    }
    if (result.hmin_effective > result.hmax_effective) {
      throw std::runtime_error(
          "computed hmin is larger than hmax. Increase hmax_fraction, "
          "reduce hmin_edge_scale, or relax cad_edge_tol_fraction.");
    }
  }

  const auto fragmentation = apply_continuity_fragmentation(
      cad, use_occ, result.mesh_size, original_face_tags);
  result.continuity_break_count = fragmentation.continuity_break_count;
  result.continuity_helper_curve_count = fragmentation.helper_curve_count;
  result.fragmented_face_count = fragmentation.fragmented_face_count;

  configure_meshing(result.mesh_size, result.hmin_effective,
                    result.hmax_effective, options);
  gmsh::model::mesh::generate(2);

  result.cad_area = occ.ok && std::isfinite(occ.surface_area) &&
                            occ.surface_area > 0.0
                        ? occ.surface_area
                        : get_occ_area();
  gmsh::model::getEntities(faces, 2);
  const auto curve_faces =
      build_curve_face_map_from_mesh(fragmentation.current_face_to_original_face);
  std::unordered_map<int, int> live_faces;
  for (const auto &face : faces) {
    if (fragmentation.current_face_to_original_face.count(face.second) != 1)
      throw std::runtime_error("[native_identity] Missing source ancestry for current face " + std::to_string(face.second));
    live_faces.emplace(face.second, face.second);
  }
  const auto current_curve_faces = build_curve_face_map_from_mesh(live_faces);
  const auto node_map = get_node_map(
      result, current_curve_faces, fragmentation.feature_curves, cad, use_occ);
  const auto ref_tris = refined_reference_triangles(options.refinement_level);
  const auto edge_nodes = legendre_gauss_lobatto_nodes(options.edge_gll_order);
  const auto edge_weights = barycentric_weights(edge_nodes);
  const auto edge_quad_weights =
      legendre_gauss_lobatto_quadrature_weights(edge_nodes);
  auto edge_helper = load_occ_edge_samples(
      cad, use_occ, edge_nodes,
      curve_faces, node_map, fragmentation.feature_curves);
  EdgeSampleMap edge_samples;
  if (edge_helper.ok) {
    edge_samples = std::move(edge_helper.edge_samples);
    result.edge_newton_success = edge_helper.newton_success;
    result.edge_curve_projection_success =
        edge_helper.curve_projection_success;
    result.edge_fallback_points = edge_helper.fallback_points;
  } else {
    edge_samples = build_edge_samples(node_map, edge_nodes, current_curve_faces,
                                      fragmentation.feature_curves, result);
  }

  auto feature_helper = load_occ_feature_samples(
      cad, use_occ, edge_nodes,
      node_map, fragmentation.feature_curves);
  if (feature_helper.ok) {
    result.feature_edge_fallback_points = feature_helper.fallback_points;
    for (auto &entry : feature_helper.edge_samples) {
      edge_samples[entry.first] = std::move(entry.second);
    }
  }
  record_edge_panels(edge_nodes, edge_weights, edge_quad_weights, edge_samples,
                     curve_faces, fragmentation.feature_curves, result);

  std::vector<Vec3> candidate_xyz;
  std::vector<int> candidate_face_tags;
  std::vector<int> candidate_current_face_tags;
  int coarse_tri_index = 0;

  for (const auto &face : faces) {
    const int current_face_tag = face.second;
    const auto original_face_it =
        fragmentation.current_face_to_original_face.find(current_face_tag);
    if (original_face_it == fragmentation.current_face_to_original_face.end()) {
      continue;
    }
    const int face_tag = original_face_it->second;
    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;
    gmsh::model::mesh::getElements(element_types, element_tags,
                                   element_node_tags, 2, current_face_tag);

    for (std::size_t block = 0; block < element_types.size(); ++block) {
      std::string element_name;
      int dim = 0;
      int element_order = 0;
      int num_nodes = 0;
      int num_primary_nodes = 0;
      std::vector<double> local_coords;
      gmsh::model::mesh::getElementProperties(
          element_types[block], element_name, dim, element_order, num_nodes,
          local_coords, num_primary_nodes);

      if (dim != 2 || num_primary_nodes < 3 ||
          element_name.find("Triangle") == std::string::npos) {
        continue;
      }

      const auto &tags = element_node_tags[block];
      const std::size_t num_elements =
          tags.size() / static_cast<std::size_t>(num_nodes);

      for (std::size_t elem = 0; elem < num_elements; ++elem) {
        const auto tag0 = tags[elem * static_cast<std::size_t>(num_nodes)];
        const auto tag1 =
            tags[elem * static_cast<std::size_t>(num_nodes) + 1];
        const auto tag2 =
            tags[elem * static_cast<std::size_t>(num_nodes) + 2];

        const auto p0_it = node_map.find(tag0);
        const auto p1_it = node_map.find(tag1);
        const auto p2_it = node_map.find(tag2);
        if (p0_it == node_map.end() || p1_it == node_map.end() ||
            p2_it == node_map.end()) {
          continue;
        }

        const Vec3 p0 = p0_it->second.xyz;
        const Vec3 p1 = p1_it->second.xyz;
        const Vec3 p2 = p2_it->second.xyz;
        const auto edge12_it = edge_samples.find(EdgeKey{tag0, tag1});
        const auto edge23_it = edge_samples.find(EdgeKey{tag1, tag2});
        const auto edge31_it = edge_samples.find(EdgeKey{tag2, tag0});
        const std::vector<Vec3> *edge12 =
            edge12_it == edge_samples.end() ? nullptr : &edge12_it->second;
        const std::vector<Vec3> *edge23 =
            edge23_it == edge_samples.end() ? nullptr : &edge23_it->second;
        const std::vector<Vec3> *edge31 =
            edge31_it == edge_samples.end() ? nullptr : &edge31_it->second;
        ++coarse_tri_index;
        result.scaffold_triangles.push_back(p0_it->second.local_index);
        result.scaffold_triangles.push_back(p1_it->second.local_index);
        result.scaffold_triangles.push_back(p2_it->second.local_index);
        result.scaffold_triangle_node_tags.push_back(static_cast<int>(tag0));
        result.scaffold_triangle_node_tags.push_back(static_cast<int>(tag1));
        result.scaffold_triangle_node_tags.push_back(static_cast<int>(tag2));
        result.scaffold_face_tags.push_back(face_tag);
        const TriangleQuality quality = compute_triangle_quality(p0, p1, p2);
        result.scaffold_min_angle_deg.push_back(quality.min_angle_deg);
        result.scaffold_edge_ratio.push_back(quality.edge_ratio);
        result.scaffold_triangle_area.push_back(quality.area);
        if (coarse_tri_index == 1 ||
            quality.min_angle_deg < result.min_scaffold_angle_deg) {
          result.min_scaffold_angle_deg = quality.min_angle_deg;
          result.worst_quality_triangle = coarse_tri_index;
        }
        result.max_scaffold_edge_ratio =
            std::max(result.max_scaffold_edge_ratio, quality.edge_ratio);
        if (quality.min_angle_deg < options.min_angle_deg) {
          ++result.bad_min_angle_count;
        }
        if (quality.edge_ratio > options.max_edge_ratio) {
          ++result.bad_edge_ratio_count;
        }

        for (const auto &subtri : ref_tris) {
          result.tri_face_tags.push_back(face_tag);
          result.parent_coarse_tri.push_back(coarse_tri_index);

          for (int node = 0; node < sample_nodes.nodes_per_patch; ++node) {
            const double u =
                sample_nodes.uv[2 * static_cast<std::size_t>(node)];
            const double v =
                sample_nodes.uv[2 * static_cast<std::size_t>(node) + 1];
            // fmm3dbie/koorn RV nodes use x(u,v) = (1-u-v) P1 + u P2 + v P3.
            const Vec3 local_bary{1.0 - u - v, u, v};
            const Vec3 parent_bary =
                map_local_bary_to_parent(local_bary, subtri);
            const Vec3 candidate_point = evaluate_blended_triangle_point(
                parent_bary, p0, p1, p2, edge_nodes, edge_weights, edge12,
                edge23, edge31);
            candidate_xyz.push_back(candidate_point);
            candidate_face_tags.push_back(face_tag);
            candidate_current_face_tags.push_back(current_face_tag);
          }
        }
      }
    }
  }

  result.nscaffold_triangles =
      static_cast<int>(result.scaffold_face_tags.size());
  if (options.enforce_quality &&
      (result.bad_min_angle_count > 0 || result.bad_edge_ratio_count > 0)) {
    std::ostringstream msg;
    msg << "Scaffold mesh quality check failed: "
        << result.bad_min_angle_count << " triangles below min_angle_deg="
        << options.min_angle_deg << ", " << result.bad_edge_ratio_count
        << " triangles above max_edge_ratio=" << options.max_edge_ratio
        << ". Worst triangle id " << result.worst_quality_triangle
        << " has min angle " << result.min_scaffold_angle_deg
        << " deg; max edge ratio is " << result.max_scaffold_edge_ratio
        << ".";
    throw std::runtime_error(msg.str());
  }

  double sum_projection_distance = 0.0;
  int projection_count = 0;
  auto projected = load_occ_projected_nodes(
      cad, use_occ, candidate_xyz,
      candidate_face_tags);
  if (projected.ok) {
    result.projection_failures = projected.fallback_points;
    for (std::size_t i = 0; i < projected.xyz.size(); ++i) {
      const double dist = projected.dist[i];
      result.proj_dist.push_back(dist);
      result.max_projection_distance =
          std::max(result.max_projection_distance, dist);
      sum_projection_distance += dist;
      ++projection_count;
      append_xyz(result.xyz, projected.xyz[i]);
    }
  } else {
    for (std::size_t i = 0; i < candidate_xyz.size(); ++i) {
      Vec3 projected_xyz = candidate_xyz[i];
      Vec2 surf_uv{};
      bool projected_ok = closest_point_on_surface(
          candidate_current_face_tags[i], candidate_xyz[i], projected_xyz, surf_uv);
      if (!projected_ok) {
        projected_ok = surface_uv_from_xyz(candidate_current_face_tags[i],
                                           candidate_xyz[i], surf_uv) &&
                       xyz_from_surface_uv(candidate_current_face_tags[i], surf_uv,
                                           projected_xyz);
      }
      if (!projected_ok) {
        ++result.projection_failures;
        projected_xyz = candidate_xyz[i];
      }

      const double dist = norm(subtract(projected_xyz, candidate_xyz[i]));
      result.proj_dist.push_back(dist);
      result.max_projection_distance =
          std::max(result.max_projection_distance, dist);
      sum_projection_distance += dist;
      ++projection_count;
      append_xyz(result.xyz, projected_xyz);
    }
  }

  result.npatches = static_cast<int>(result.tri_face_tags.size());
  result.mean_projection_distance =
      projection_count > 0 ? sum_projection_distance / projection_count : 0.0;
  return result;
}

} // namespace stepmesher

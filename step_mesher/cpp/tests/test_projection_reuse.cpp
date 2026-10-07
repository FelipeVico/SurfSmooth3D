// Differential regression and bounded microbenchmark for projector reuse.
// The oracle below deliberately retains the pre-reuse projection semantics.
// No meshing or Gmsh model construction is performed by this test.
#include "step_mesher/occ_context.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_BSplineSurface.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp.hxx>
#include <TopTools_FormatVersion.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
using stepmesher::occ::Context;
using stepmesher::occ::ProjectNode;
using stepmesher::occ::ProjectedNodes;
using stepmesher::occ::Vec2;
using stepmesher::occ::Vec3;
using Clock = std::chrono::steady_clock;

void require(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}
bool finite(const Vec3 &x) {
  return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}
Vec3 pnt_to_vec(const gp_Pnt &p) { return {p.X(), p.Y(), p.Z()}; }
Vec3 subtract(const Vec3 &a, const Vec3 &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
double norm(const Vec3 &a) {
  return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

// Frozen copy of occ_context.cpp::project_surface before projector reuse.
bool frozen_fresh_project_surface(const TopoDS_Face &face, const Vec3 &x,
                                  Vec2 &uv, Vec3 &out) {
  try {
    Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
    if (surf.IsNull()) return false;
    BRepAdaptor_Surface adapt(face, Standard_False);
    gp_Pnt p(x[0], x[1], x[2]);
    GeomAPI_ProjectPointOnSurf proj;
    const double u0 = adapt.FirstUParameter();
    const double u1 = adapt.LastUParameter();
    const double v0 = adapt.FirstVParameter();
    const double v1 = adapt.LastVParameter();
    if (std::isfinite(u0) && std::isfinite(u1) && std::isfinite(v0) &&
        std::isfinite(v1) && std::abs(u1 - u0) > 1.0e-30 &&
        std::abs(v1 - v0) > 1.0e-30 && std::abs(u0) < 1.0e50 &&
        std::abs(u1) < 1.0e50 && std::abs(v0) < 1.0e50 &&
        std::abs(v1) < 1.0e50) {
      proj.Init(p, surf, u0, u1, v0, v1);
    } else {
      proj.Init(p, surf);
    }
    if (!proj.IsDone() || proj.NbPoints() <= 0) return false;
    Standard_Real u = 0.0;
    Standard_Real v = 0.0;
    proj.LowerDistanceParameters(u, v);
    uv = {u, v};
    out = pnt_to_vec(proj.NearestPoint());
    return finite(out);
  } catch (...) {
    return false;
  }
}

ProjectedNodes frozen_fresh_project_nodes(const std::vector<TopoDS_Face> &faces,
                                         const std::vector<ProjectNode> &nodes) {
  ProjectedNodes result;
  for (const ProjectNode &node : nodes) {
    require(node.face_tag > 0 &&
                static_cast<std::size_t>(node.face_tag) <= faces.size(),
            "invalid oracle face tag");
    Vec3 projected = node.x;
    Vec2 uv{};
    const bool ok = frozen_fresh_project_surface(
        faces[static_cast<std::size_t>(node.face_tag - 1)], node.x, uv, projected);
    if (!ok) {
      ++result.n_fallback;
      uv = {0.0, 0.0};
      projected = node.x;
    }
    const double dist = ok ? norm(subtract(projected, node.x)) : 0.0;
    result.nodes.push_back({projected, dist, uv});
  }
  return result;
}

std::vector<TopoDS_Face> snapshot_faces(Context &context, const std::string &step_file) {
  const std::string snapshot = context.adapter_brep_snapshot();
  std::istringstream input(snapshot);
  BRep_Builder builder;
  TopoDS_Shape shape;
  BRepTools::Read(shape, input, builder);
  require(!shape.IsNull(), "could not read retained BRep snapshot");
  std::ostringstream roundtrip;
  BRepTools::Write(shape, roundtrip, Standard_False, Standard_False,
                  TopTools_FormatVersion_CURRENT);
  if (roundtrip.str() != snapshot) {
    // Reading analytic axes can renormalize already-normalized directions and
    // change their last bits. Do not certify that altered shape as the oracle.
    // A separate STEP import avoids this round trip; accept it only when its
    // complete serialization exactly matches the retained Context geometry.
    STEPControl_Reader reader;
    require(reader.ReadFile(step_file.c_str()) == IFSelect_RetDone,
            "could not import unroundtripped oracle STEP geometry");
    reader.TransferRoots();
    shape = reader.OneShape();
    require(!shape.IsNull(), "unroundtripped oracle STEP import is empty");
    std::ostringstream imported;
    BRepTools::Write(shape, imported, Standard_False, Standard_False,
                    TopTools_FormatVersion_CURRENT);
    require(imported.str() == snapshot,
            "neither BRep roundtrip nor independent STEP import exactly "
            "matches retained geometry; oracle unverified");
  }
  TopTools_IndexedMapOfShape face_map;
  TopExp::MapShapes(shape, TopAbs_FACE, face_map);
  require(static_cast<std::size_t>(face_map.Extent()) == context.face_count(),
          "retained BRep face count changed");
  std::vector<TopoDS_Face> faces;
  std::unordered_map<int, int> mapping;
  for (int i = 1; i <= face_map.Extent(); ++i) {
    faces.push_back(TopoDS::Face(face_map.FindKey(i)));
    mapping.emplace(i, i - 1);
  }
  context.bind_adapter_face_indices(mapping);
  return faces;
}

std::string json_string(const std::string &value) {
  std::ostringstream out;
  out << '"';
  for (unsigned char c : value) {
    if (c == '"' || c == '\\') out << '\\' << c;
    else if (c == '\n') out << "\\n";
    else if (c == '\r') out << "\\r";
    else if (c == '\t') out << "\\t";
    else if (c < 32) out << "\\u" << std::hex << std::setw(4)
                         << std::setfill('0') << static_cast<int>(c) << std::dec;
    else out << c;
  }
  out << '"';
  return out.str();
}

std::string digest(const ProjectedNodes &result) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  auto bytes = [&](const void *ptr, std::size_t count) {
    const auto *p = static_cast<const unsigned char *>(ptr);
    for (std::size_t i = 0; i < count; ++i) {
      hash ^= p[i];
      hash *= UINT64_C(1099511628211);
    }
  };
  bytes(&result.n_fallback, sizeof(result.n_fallback));
  for (const auto &node : result.nodes) {
    for (double value : node.x) bytes(&value, sizeof(value));
    bytes(&node.distance, sizeof(node.distance));
    for (double value : node.uv) bytes(&value, sizeof(value));
  }
  std::ostringstream out;
  out << std::hex << std::setw(16) << std::setfill('0') << hash;
  return out.str();
}

void compare(const ProjectedNodes &oracle, const ProjectedNodes &candidate,
             const std::string &label) {
  require(oracle.nodes.size() == candidate.nodes.size(), label + ": node count changed");
  require(oracle.n_fallback == candidate.n_fallback, label + ": fallback count changed");
  for (std::size_t i = 0; i < oracle.nodes.size(); ++i) {
    const auto &a = oracle.nodes[i];
    const auto &b = candidate.nodes[i];
    auto same_double = [](double first, double second) {
      return std::memcmp(&first, &second, sizeof(double)) == 0;
    };
    bool same = same_double(a.distance, b.distance);
    for (std::size_t component = 0; component < a.x.size(); ++component)
      same = same && same_double(a.x[component], b.x[component]);
    for (std::size_t component = 0; component < a.uv.size(); ++component)
      same = same && same_double(a.uv[component], b.uv[component]);
    if (!same) {
      std::ostringstream error;
      error << std::setprecision(17) << label << ": exact mismatch at node " << i;
      error << "; oracle=";
      for (double value : a.x) error << value << ',';
      error << a.uv[0] << ',' << a.uv[1] << ',' << a.distance << "; candidate=";
      for (double value : b.x) error << value << ',';
      error << b.uv[0] << ',' << b.uv[1] << ',' << b.distance;
      error << "; oracle digest=" << digest(oracle)
            << "; candidate digest=" << digest(candidate);
      throw std::runtime_error(error.str());
    }
  }
}

struct ReplayInputs {
  std::string step_file;
  std::size_t source_count = 0, face_count = 0;
  std::vector<ProjectNode> nodes;
};

ReplayInputs read_replay(const std::string &request_path, const std::string &map_path,
                         std::size_t cap, std::size_t face_count) {
  std::ifstream mapping_file(map_path);
  require(static_cast<bool>(mapping_file), "could not open replay canonical face map");
  std::unordered_map<int, int> mapping;
  int tag = 0, index = 0;
  while (mapping_file >> tag >> index) {
    require(index >= 0 && static_cast<std::size_t>(index) < face_count,
            "replay map contains invalid canonical face index");
    require(mapping.emplace(tag, index).second, "replay map contains duplicate tag");
  }
  require(mapping_file.eof() && !mapping.empty(), "malformed or empty replay face map");
  std::ifstream request(request_path);
  require(static_cast<bool>(request), "could not open replay projection request");
  std::string marker;
  std::getline(request, marker);
  require(marker == "STEP_MESHER_OCC_PROJECT_V1", "bad replay projection request marker");
  ReplayInputs result;
  std::getline(request, result.step_file);
  std::size_t descriptor_count = 0;
  request >> descriptor_count;
  for (std::size_t i = 0; i < descriptor_count; ++i) {
    int descriptor_tag = 0;
    double ignored = 0;
    request >> descriptor_tag;
    for (int field = 0; field < 7; ++field) request >> ignored;
  }
  request >> result.source_count;
  std::vector<ProjectNode> source(result.source_count);
  std::map<int, std::vector<std::size_t>> groups;
  for (std::size_t i = 0; i < source.size(); ++i) {
    auto &node = source[i];
    request >> node.face_tag >> node.x[0] >> node.x[1] >> node.x[2];
    require(static_cast<bool>(request), "truncated replay projection request");
    const auto association = mapping.find(node.face_tag);
    require(association != mapping.end(), "replay point has no canonical face association");
    node.face_tag = association->second + 1;
    groups[node.face_tag].push_back(i);
  }
  require(static_cast<bool>(request) && !source.empty(), "empty replay projection request");
  result.face_count = groups.size();
  require(cap >= groups.size(), "replay cap must permit at least one point per represented face");
  if (source.size() <= cap) {
    result.nodes = std::move(source);
    return result;
  }
  // Allocate slots evenly over faces, then evenly over each face's source order.
  // Sort selected source indices to preserve the original candidate traversal.
  std::map<int, std::size_t> allocation;
  std::size_t allocated = 0;
  while (allocated < cap) {
    for (const auto &group : groups) {
      auto &count = allocation[group.first];
      if (count < group.second.size()) { ++count; ++allocated; }
      if (allocated == cap) break;
    }
  }
  std::vector<std::size_t> selected;
  selected.reserve(cap);
  for (const auto &group : groups) {
    const std::size_t count = allocation.at(group.first);
    for (std::size_t i = 0; i < count; ++i) {
      const std::size_t position = count == 1 ? group.second.size() / 2 :
          i * (group.second.size() - 1) / (count - 1);
      selected.push_back(group.second[position]);
    }
  }
  std::sort(selected.begin(), selected.end());
  result.nodes.reserve(selected.size());
  for (std::size_t position : selected) result.nodes.push_back(source[position]);
  return result;
}

std::pair<double, double> narrowest_span(const Handle(Geom_BSplineSurface) &s,
                                       bool u_direction, double lower, double upper) {
  std::pair<double, double> best{lower, upper};
  if (s.IsNull()) return best;
  const int count = u_direction ? s->NbUKnots() : s->NbVKnots();
  for (int i = 1; i < count; ++i) {
    const double a = std::max(lower, u_direction ? s->UKnot(i) : s->VKnot(i));
    const double b = std::min(upper, u_direction ? s->UKnot(i + 1) : s->VKnot(i + 1));
    if (b > a && b - a < best.second - best.first) best = {a, b};
  }
  return best;
}

std::size_t parameter_count(const Handle(Geom_BSplineSurface) &s, bool u_direction,
                            double lower, double upper) {
  const int degree = u_direction ? s->UDegree() : s->VDegree();
  const int knots = u_direction ? s->NbUKnots() : s->NbVKnots();
  double previous = lower;
  std::size_t count = 1;
  for (int i = 1; i < knots; ++i) {
    const double a = u_direction ? s->UKnot(i) : s->VKnot(i);
    const double b = u_direction ? s->UKnot(i + 1) : s->VKnot(i + 1);
    if (a >= upper - Precision::PConfusion()) break;
    if (b < lower + Precision::PConfusion()) continue;
    const double step = (b - a) / std::max(degree, 2);
    for (int k = 1; k <= degree; ++k) {
      const double parameter = a + k * step;
      if (parameter > upper - Precision::PConfusion()) break;
      if (parameter > previous + Precision::PConfusion()) {
        ++count;
        previous = parameter;
      }
    }
  }
  return std::max<std::size_t>(44, count + 1);
}

std::size_t grid_estimate(const TopoDS_Face &face) {
  BRepAdaptor_Surface adapt(face, Standard_False);
  if (adapt.GetType() != GeomAbs_BSplineSurface) return 0;
  const auto spline = adapt.BSpline();
  return parameter_count(spline, true, adapt.FirstUParameter(), adapt.LastUParameter()) *
         parameter_count(spline, false, adapt.FirstVParameter(), adapt.LastVParameter());
}

const char *surface_type(const TopoDS_Face &face) {
  BRepAdaptor_Surface adapt(face, Standard_False);
  switch (adapt.GetType()) {
    case GeomAbs_Plane: return "plane";
    case GeomAbs_Cylinder: return "cylinder";
    case GeomAbs_Cone: return "cone";
    case GeomAbs_Sphere: return "sphere";
    case GeomAbs_Torus: return "torus";
    case GeomAbs_BSplineSurface: return "bspline";
    case GeomAbs_BezierSurface: return "bezier";
    case GeomAbs_SurfaceOfRevolution: return "revolution";
    case GeomAbs_SurfaceOfExtrusion: return "extrusion";
    case GeomAbs_OffsetSurface: return "offset";
    default: return "other";
  }
}

std::vector<ProjectNode> sample_face(const TopoDS_Face &face, int tag,
                                     std::size_t count, double scale) {
  BRepAdaptor_Surface adapt(face, Standard_False);
  double u0 = adapt.FirstUParameter(), u1 = adapt.LastUParameter();
  double v0 = adapt.FirstVParameter(), v1 = adapt.LastVParameter();
  if (!std::isfinite(u0) || !std::isfinite(u1) || !std::isfinite(v0) ||
      !std::isfinite(v1) || std::max({std::abs(u0), std::abs(u1),
                                    std::abs(v0), std::abs(v1)}) >= 1.0e50) {
    BRepTools::UVBounds(face, u0, u1, v0, v1);
  }
  require(std::isfinite(u0) && std::isfinite(u1) && std::isfinite(v0) &&
              std::isfinite(v1) && u1 > u0 && v1 > v0,
          "cannot sample face domain");
  Handle(Geom_BSplineSurface) spline;
  if (adapt.GetType() == GeomAbs_BSplineSurface) spline = adapt.BSpline();
  const auto uspan = narrowest_span(spline, true, u0, u1);
  const auto vspan = narrowest_span(spline, false, v0, v1);
  const double unarrow = 0.5 * (uspan.first + uspan.second);
  const double vnarrow = 0.5 * (vspan.first + vspan.second);
  const auto u = [&](double fraction) { return u0 + fraction * (u1 - u0); };
  const auto v = [&](double fraction) { return v0 + fraction * (v1 - v0); };
  const std::vector<Vec2> parameters = {
      {u(0.173), v(0.337)}, {u(0.173), v(0.337)},
      {u0, v(0.413)}, {u1, v(0.413)},
      {u(0.271), v0}, {u(0.271), v1},
      {unarrow, v(0.417)}, {u(0.311), vnarrow},
      {unarrow, vnarrow}, {u0, v0},
      {u(0.739), v(0.683)}, {u(0.5), v(0.5)}};
  std::vector<ProjectNode> out;
  out.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    // Subsequent cycles also probe just inside the seams and alter the offset.
    Vec2 uv = parameters[i % parameters.size()];
    if (i >= parameters.size()) {
      const double delta = 1.0e-7 * static_cast<double>(i / parameters.size());
      uv[0] = std::clamp(uv[0] + delta * (u1 - u0), u0, u1);
      uv[1] = std::clamp(uv[1] - delta * (v1 - v0), v0, v1);
    }
    gp_Pnt point;
    gp_Vec du, dv;
    adapt.D1(uv[0], uv[1], point, du, dv);
    gp_Vec normal = du.Crossed(dv);
    if (normal.SquareMagnitude() > 1.0e-30) normal.Normalize();
    else normal = gp_Vec(1.0, -0.5, 0.25).Normalized();
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    const double offset = (i % 3 == 2) ? 0.0 : sign * scale * 1.0e-4;
    point.Translate(normal * offset);
    const Vec3 xyz = pnt_to_vec(point);
    require(finite(xyz), "generated nonfinite projection input");
    out.push_back({tag, xyz});
  }
  return out;
}

template <class Operation>
std::pair<ProjectedNodes, double> timed(Operation operation) {
  const auto start = Clock::now();
  auto result = operation();
  const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
  return {std::move(result), seconds};
}

void check_contracts(Context &context, const std::string &fixture,
                     const std::vector<TopoDS_Face> &faces,
                     std::size_t primary, double scale) {
  const auto empty = context.project_nodes({});
  require(empty.nodes.empty() && empty.n_fallback == 0,
          "empty projection batch changed behavior");
  const auto valid = sample_face(faces[primary], static_cast<int>(primary + 1), 2, scale);
  const auto expected = frozen_fresh_project_nodes(faces, valid);
  compare(expected, context.project_nodes(valid), "valid batch before contract rejection");
  auto mixed = valid;
  ProjectNode unknown = valid.front();
  unknown.face_tag = static_cast<int>(faces.size() + 1);
  mixed.insert(mixed.begin() + 1, unknown);
  bool rejected = false;
  try {
    context.project_nodes(mixed);
  } catch (const stepmesher::occ::NativeContractError &error) {
    rejected = std::string(error.what()).rfind("[native_identity]", 0) == 0;
  }
  require(rejected, "unknown face must raise NativeContractError before projection");
  compare(expected, context.project_nodes(valid), "valid batch after contract rejection");

  Context independent(fixture);
  const auto independent_faces = snapshot_faces(independent, fixture);
  require(independent_faces.size() == faces.size(), "independent context face count changed");
  compare(expected, independent.project_nodes(valid), "independent context");
  compare(expected, context.project_nodes(valid), "original context after independent context");
  std::cout << "{\"fixture\":" << json_string(fixture)
            << ",\"case\":\"contracts\",\"empty_batch\":true,"
               "\"unknown_face_rejected\":true,\"independent_contexts\":true,"
               "\"exact\":true}\n" << std::flush;
}

void check_failure_recovery(Context &context, const std::string &fixture,
                            const std::vector<TopoDS_Face> &faces,
                            std::size_t preferred, double scale) {
  std::size_t face_index = preferred;
  auto spline = [&](std::size_t index) {
    BRepAdaptor_Surface adapt(faces[index], Standard_False);
    return adapt.GetType() == GeomAbs_BSplineSurface ||
           adapt.GetType() == GeomAbs_BezierSurface;
  };
  if (!spline(face_index)) {
    face_index = 0;
    while (face_index < faces.size() && !spline(face_index)) ++face_index;
  }
  require(face_index < faces.size(), "failure probe needs a BSpline or Bezier face");
  const auto valid = sample_face(faces[face_index], static_cast<int>(face_index + 1), 1, scale);
  auto poisoned = valid.front();
  poisoned.x[0] = std::numeric_limits<double>::quiet_NaN();
  const std::vector<ProjectNode> sequence{valid.front(), poisoned, valid.front()};
  auto phase = [&](const char *stage) {
    std::cout << "{\"fixture\":" << json_string(fixture)
              << ",\"case\":\"failure probe phase\",\"stage\":"
              << json_string(stage) << "}\n" << std::flush;
  };
  phase("fresh oracle");
  const auto expected = frozen_fresh_project_nodes(faces, sequence);
  require(expected.n_fallback == 1, "nonfinite oracle point did not produce one fallback");
  phase("candidate");
  const auto actual = context.project_nodes(sequence);
  compare(expected, actual, "valid/nonfinite/valid failure recovery");
  const auto repeated = context.project_nodes(sequence);
  compare(expected, repeated, "repeated valid/nonfinite/valid failure recovery");
  compare(frozen_fresh_project_nodes(faces, valid), context.project_nodes(valid),
          "valid point after failed batch");
  std::cout << "{\"fixture\":" << json_string(fixture)
            << ",\"case\":\"failure recovery\",\"face\":" << face_index + 1
            << ",\"fallbacks\":" << actual.n_fallback
            << ",\"oracle_digest\":" << json_string(digest(expected))
            << ",\"candidate_digest\":" << json_string(digest(actual))
            << ",\"exact\":true}\n";
}

void report(const std::string &fixture, const std::string &label, int face,
            const char *type, std::size_t estimate, std::size_t count,
            double oracle_seconds, double candidate_seconds,
            const ProjectedNodes &oracle, const ProjectedNodes &candidate) {
  std::cout << std::setprecision(17)
            << "{\"fixture\":" << json_string(fixture)
            << ",\"case\":" << json_string(label)
            << ",\"face\":" << face << ",\"surface_type\":" << json_string(type)
            << ",\"knot_grid_points_estimate\":" << estimate
            << ",\"nodes\":" << count
            << ",\"oracle_seconds\":" << oracle_seconds
            << ",\"candidate_seconds\":" << candidate_seconds
            << ",\"speedup\":" << (candidate_seconds > 0.0 ? oracle_seconds / candidate_seconds : 0.0)
            << ",\"fallbacks\":" << candidate.n_fallback
            << ",\"oracle_digest\":" << json_string(digest(oracle))
            << ",\"candidate_digest\":" << json_string(digest(candidate))
            << ",\"exact\":true}\n";
}

std::vector<std::size_t> parse_sizes(const std::string &text) {
  std::vector<std::size_t> sizes;
  std::istringstream input(text);
  std::string word;
  while (std::getline(input, word, ',')) {
    std::size_t used = 0;
    const auto count = std::stoul(word, &used);
    require(used == word.size() && count > 0 && count <= 4096,
            "batch sizes must be integers from 1 through 4096");
    sizes.push_back(count);
  }
  require(!sizes.empty(), "no batch sizes supplied");
  return sizes;
}
} // namespace

int main(int argc, char **argv) {
  try {
    require(argc >= 2,
            "usage: test_projection_reuse STEP [--batch-sizes=1,2,16,64] "
            "[--points-per-face=12] [--repetitions=1] "
            "[--request=PATH --map=PATH --replay-cap=256] [--failure-probe]");
    const std::string fixture = argv[1];
    std::vector<std::size_t> sizes{1, 2, 16, 64};
    std::size_t points_per_face = 12;
    std::size_t repetitions = 1;
    std::size_t replay_cap = 256;
    std::string request_path, map_path;
    bool failure_probe = false;
    for (int i = 2; i < argc; ++i) {
      const std::string argument = argv[i];
      if (argument.rfind("--batch-sizes=", 0) == 0)
        sizes = parse_sizes(argument.substr(14));
      else if (argument.rfind("--points-per-face=", 0) == 0)
        points_per_face = parse_sizes(argument.substr(18)).at(0);
      else if (argument.rfind("--repetitions=", 0) == 0)
        repetitions = parse_sizes(argument.substr(14)).at(0);
      else if (argument.rfind("--request=", 0) == 0)
        request_path = argument.substr(10);
      else if (argument.rfind("--map=", 0) == 0)
        map_path = argument.substr(6);
      else if (argument.rfind("--replay-cap=", 0) == 0)
        replay_cap = parse_sizes(argument.substr(13)).at(0);
      else if (argument == "--failure-probe")
        failure_probe = true;
      else throw std::runtime_error("unknown argument: " + argument);
    }
    Context context(fixture);
    const std::string original_snapshot = context.adapter_brep_snapshot();
    const auto faces = snapshot_faces(context, fixture);
    require(!faces.empty(), "fixture has no faces");
    require(request_path.empty() == map_path.empty(),
            "replay requires both --request and --map");
    const double scale = std::max(1.0, context.measurements(1.0e-6).diag);
    std::size_t primary = 0, largest_grid = 0;
    std::vector<ProjectNode> all_faces;
    for (std::size_t i = 0; i < faces.size(); ++i) {
      const auto estimate = grid_estimate(faces[i]);
      if (estimate > largest_grid) {
        primary = i;
        largest_grid = estimate;
      }
      const auto inputs = sample_face(faces[i], static_cast<int>(i + 1), points_per_face, scale);
      all_faces.insert(all_faces.end(), inputs.begin(), inputs.end());
    }
    check_contracts(context, fixture, faces, primary, scale);
    if (failure_probe) {
      require(request_path.empty(), "--failure-probe cannot be combined with replay");
      check_failure_recovery(context, fixture, faces, primary, scale);
      require(context.adapter_brep_snapshot() == original_snapshot,
              "failure probe changed retained geometry");
      return 0;
    }
    const bool replay = !request_path.empty();
    if (replay) {
      auto inputs = read_replay(request_path, map_path, replay_cap, faces.size());
      all_faces = std::move(inputs.nodes);
      std::cout << "{\"fixture\":" << json_string(fixture)
                << ",\"case\":\"replay selection\",\"request\":" << json_string(request_path)
                << ",\"request_step_file\":" << json_string(inputs.step_file)
                << ",\"source_nodes\":" << inputs.source_count
                << ",\"selected_nodes\":" << all_faces.size()
                << ",\"represented_faces\":" << inputs.face_count
                << ",\"selection\":\"evenly over faces then source order\"}\n";
    }
    const std::string prefix = replay ? "replay" : "all-faces";

    // Mixed-face A/B/A checks catch state leaking between faces or batches.
    auto a = timed([&] { return frozen_fresh_project_nodes(faces, all_faces); });
    auto candidate_a = timed([&] { return context.project_nodes(all_faces); });
    compare(a.first, candidate_a.first, prefix + " A");
    report(fixture, prefix + " A", 0, "mixed", 0, all_faces.size(),
           a.second, candidate_a.second, a.first, candidate_a.first);
    std::vector<ProjectNode> shuffled = all_faces;
    std::mt19937 random(562968);
    std::shuffle(shuffled.begin(), shuffled.end(), random);
    auto b = timed([&] { return frozen_fresh_project_nodes(faces, shuffled); });
    auto candidate_b = timed([&] { return context.project_nodes(shuffled); });
    compare(b.first, candidate_b.first, prefix + " shuffled B");
    report(fixture, prefix + " shuffled B", 0, "mixed", 0, shuffled.size(),
           b.second, candidate_b.second, b.first, candidate_b.first);
    auto repeated_a = timed([&] { return context.project_nodes(all_faces); });
    compare(a.first, repeated_a.first, prefix + " repeated A");
    report(fixture, prefix + " repeated A", 0, "mixed", 0, all_faces.size(),
           a.second, repeated_a.second, a.first, repeated_a.first);

    if (!replay) for (std::size_t count : sizes) {
      const auto inputs = sample_face(faces[primary], static_cast<int>(primary + 1), count, scale);
      for (std::size_t repeat = 0; repeat < repetitions; ++repeat) {
        auto oracle = timed([&] { return frozen_fresh_project_nodes(faces, inputs); });
        auto candidate = timed([&] { return context.project_nodes(inputs); });
        const std::string label = "primary batch " + std::to_string(count) +
                                  " repetition " + std::to_string(repeat + 1);
        compare(oracle.first, candidate.first, label);
        report(fixture, label, static_cast<int>(primary + 1), surface_type(faces[primary]),
               largest_grid, count, oracle.second, candidate.second, oracle.first, candidate.first);
      }
    }
    require(context.adapter_brep_snapshot() == original_snapshot,
            "projection changed retained geometry");
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "test_projection_reuse: " << error.what() << '\n';
    return 1;
  } catch (...) {
    std::cerr << "test_projection_reuse: unknown OCCT failure\n";
    return 1;
  }
}

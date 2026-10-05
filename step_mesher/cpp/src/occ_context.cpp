#include "step_mesher/occ_context.hpp"
#include <gmsh.h>
#include <Interface_Static.hxx>
#include <STEPControl_Controller.hxx>
#include <Standard_Version.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepTools.hxx>
#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <Bnd_Box.hxx>
#include <Extrema_ExtPC.hxx>
#include <GeomAbs_Shape.hxx>
#include <GProp_GProps.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_Handle.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>


namespace stepmesher::occ {
namespace {
using FaceDesc = FaceDescriptor;
using Request = EdgeRequest;

struct OccFace {
  TopoDS_Face face;
};

struct OccEdge {
  TopoDS_Edge edge;
  bool degenerated = false;
  std::vector<int> face_ids;
};


struct SurfaceEval {
  Vec3 x{};
  Vec3 xu{};
  Vec3 xv{};
  Vec3 xuu{};
  Vec3 xvv{};
  Vec3 xuv{};
};


double edge_length(const TopoDS_Edge &edge) {
  GProp_GProps props;
  BRepGProp::LinearProperties(edge, props);
  return props.Mass();
}

double norm(const Vec3 &a) {
  return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

double dot(const Vec3 &a, const Vec3 &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 cross(const Vec3 &a, const Vec3 &b) {
  return {a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

Vec3 add(const Vec3 &a, const Vec3 &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vec3 subtract(const Vec3 &a, const Vec3 &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec3 scale(const Vec3 &a, double s) {
  return {s * a[0], s * a[1], s * a[2]};
}

Vec3 pnt_to_vec(const gp_Pnt &p) {
  return {p.X(), p.Y(), p.Z()};
}

Vec3 vec_to_vec(const gp_Vec &v) {
  return {v.X(), v.Y(), v.Z()};
}

bool finite(const Vec3 &a) {
  return std::isfinite(a[0]) && std::isfinite(a[1]) && std::isfinite(a[2]);
}

FaceDesc make_desc(const std::array<double, 6> &bbox, double area, int tag) {
  FaceDesc d;
  d.tag = tag;
  d.bbox = bbox;
  d.center = {0.5 * (bbox[0] + bbox[3]), 0.5 * (bbox[1] + bbox[4]),
              0.5 * (bbox[2] + bbox[5])};
  d.size = {bbox[3] - bbox[0], bbox[4] - bbox[1], bbox[5] - bbox[2]};
  d.diag = norm(d.size);
  d.area = area;
  return d;
}

FaceDesc face_desc_occ(const TopoDS_Face &face, int tag) {
  Bnd_Box box;
  BRepBndLib::Add(face, box);
  Standard_Real x0 = 0.0;
  Standard_Real y0 = 0.0;
  Standard_Real z0 = 0.0;
  Standard_Real x1 = 0.0;
  Standard_Real y1 = 0.0;
  Standard_Real z1 = 0.0;
  box.Get(x0, y0, z0, x1, y1, z1);

  GProp_GProps props;
  BRepGProp::SurfaceProperties(face, props, 1.0e-9);
  return make_desc({x0, y0, z0, x1, y1, z1}, props.Mass(), tag);
}

double bbox_diag(const TopoDS_Shape &shape) {
  Bnd_Box box;
  BRepBndLib::AddOptimal(shape, box, Standard_False, Standard_False);

  Standard_Real x0 = 0.0;
  Standard_Real y0 = 0.0;
  Standard_Real z0 = 0.0;
  Standard_Real x1 = 0.0;
  Standard_Real y1 = 0.0;
  Standard_Real z1 = 0.0;
  box.Get(x0, y0, z0, x1, y1, z1);

  return norm(Vec3{x1 - x0, y1 - y0, z1 - z0});
}


// Freeze the defaults of the old, separately launched STEP reader. Its option
// fields were never applied; this scope prevents Gmsh/prior callers from leaking
// process-global reader changes into the now in-process reader.
class ReaderSettings {
public:
  ReaderSettings() {
    STEPControl_Controller::Init();
    const std::pair<const char *, const char *> defaults[] = {
      {"read.precision.mode", "File"}, {"read.precision.val", "0.001"},
      {"read.maxprecision.mode", "Preferred"}, {"read.maxprecision.val", "1"},
      {"read.surfacecurve.mode", "Default"},
      {"read.stdsameparameter.mode", "Off"}, {"xstep.cascade.unit", "MM"},
      {"read.step.shape.repr", "All"}, {"read.step.shape.relationship", "ON"},
      {"read.step.shape.aspect", "ON"}, {"read.step.product.mode", "ON"},
      {"read.step.product.context", "all"}, {"read.step.assembly.level", "all"},
      {"read.step.nonmanifold", "Off"}, {"read.step.ideas", "Off"},
      {"read.step.resource.name", "STEP"}, {"read.step.sequence", "FromSTEP"},
      {"read.encoderegularity.angle", "0.01"}, {"step.angleunit.mode", "File"},
      {"read.step.tessellated", "On"},
      {"read.step.constructivegeom.relationship", "OFF"},
      {"read.step.codepage", "UTF8"}, {"read.step.all.shapes", "Off"},
      {"read.step.root.transformation", "ON"}
    };
    for (const auto &setting : defaults) {
      const char *old = Interface_Static::CVal(setting.first);
      previous_.emplace_back(setting.first, old ? old : "");
      if (!Interface_Static::SetCVal(setting.first, setting.second)) {
        restore();
        throw NativeContractError(std::string("[native_runtime] unsupported reader setting ") + setting.first);
      }
    }
    // These statics are absent in a fresh STEPControl-only helper; IVal then
    // returns zero. Other CAD clients can register them in this shared process.
    for (const char *name : {"read.iges.bspline.continuity", "read.stepcaf.subshapes.name"}) {
      const char *old = Interface_Static::CVal(name);
      if (old && *old) {
        previous_.emplace_back(name, old);
        if (!Interface_Static::SetIVal(name, 0)) {
          restore();
          throw NativeContractError(std::string("[native_runtime] could not freeze optional reader setting ") + name);
        }
      }
    }
  }
  ~ReaderSettings() { restore(); }
private:
  void restore() noexcept {
    for (auto it = previous_.rbegin(); it != previous_.rend(); ++it)
      Interface_Static::SetCVal(it->first.c_str(), it->second.c_str());
    previous_.clear();
  }
  std::vector<std::pair<std::string, std::string>> previous_;
};

TopoDS_Shape read_step(const std::string &step_file) {
  ReaderSettings settings;
  STEPControl_Reader reader;
  const IFSelect_ReturnStatus status = reader.ReadFile(step_file.c_str());
  if (status != IFSelect_RetDone) {
    throw std::runtime_error("could not read STEP file");
  }
  reader.TransferRoots();
  TopoDS_Shape shape = reader.OneShape();
  if (shape.IsNull()) {
    throw std::runtime_error("STEP file produced an empty shape");
  }
  return shape;
}


std::vector<OccFace> collect_faces(const TopoDS_Shape &shape) {
  TopTools_IndexedMapOfShape face_map;
  TopExp::MapShapes(shape, TopAbs_FACE, face_map);
  std::vector<OccFace> faces;
  faces.reserve(static_cast<std::size_t>(face_map.Extent()));
  for (int fi = 1; fi <= face_map.Extent(); ++fi) {
    const TopoDS_Face face = TopoDS::Face(face_map.FindKey(fi));
    faces.push_back({face});
  }
  return faces;
}

std::vector<OccEdge> collect_edges(const TopoDS_Shape &shape) {
  TopTools_IndexedMapOfShape face_map;
  TopTools_IndexedMapOfShape edge_map;
  TopExp::MapShapes(shape, TopAbs_FACE, face_map);
  TopExp::MapShapes(shape, TopAbs_EDGE, edge_map);

  std::vector<std::vector<int>> edge_faces(edge_map.Extent() + 1);
  for (int fi = 1; fi <= face_map.Extent(); ++fi) {
    TopoDS_Face face = TopoDS::Face(face_map.FindKey(fi));
    for (TopExp_Explorer exp(face, TopAbs_EDGE); exp.More(); exp.Next()) {
      TopoDS_Edge edge = TopoDS::Edge(exp.Current());
      const int ei = edge_map.FindIndex(edge);
      if (ei > 0) {
        edge_faces[static_cast<std::size_t>(ei)].push_back(fi - 1);
      }
    }
  }

  std::vector<OccEdge> edges;
  for (int ei = 1; ei <= edge_map.Extent(); ++ei) {
    TopoDS_Edge edge = TopoDS::Edge(edge_map.FindKey(ei));
    bool deg = false;
    try {
      deg = BRep_Tool::Degenerated(edge);
    } catch (...) {
    }
    auto faces = edge_faces[static_cast<std::size_t>(ei)];
    std::sort(faces.begin(), faces.end());
    faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
    edges.push_back({edge, deg, faces});
  }
  return edges;
}

bool project_curve(const Vec3 &x, const TopoDS_Edge &edge, Vec3 &out) {
  try {
    BRepAdaptor_Curve curve(edge);
    gp_Pnt p(x[0], x[1], x[2]);
    Extrema_ExtPC ext(p, curve);
    if (!ext.IsDone() || ext.NbExt() <= 0) {
      return false;
    }
    double best = std::numeric_limits<double>::infinity();
    gp_Pnt best_p;
    for (int k = 1; k <= ext.NbExt(); ++k) {
      const double d2 = ext.SquareDistance(k);
      if (d2 < best) {
        best = d2;
        best_p = ext.Point(k).Value();
      }
    }
    out = pnt_to_vec(best_p);
    return finite(out);
  } catch (...) {
    return false;
  }
}

double point_curve_distance(const Vec3 &x, const TopoDS_Edge &edge) {
  Vec3 out{};
  if (!project_curve(x, edge, out)) {
    return std::numeric_limits<double>::infinity();
  }
  return norm(subtract(x, out));
}

int best_occ_edge_for_curve(
    int curve_tag, const Request &req, const std::vector<OccEdge> &edges,
    const std::vector<int> &preferred_faces) {
  std::vector<Vec3> samples;
  for (const auto &line : req.lines) {
    if (line.curve != curve_tag) {
      continue;
    }
    const auto a = req.nodes.find(line.a);
    const auto b = req.nodes.find(line.b);
    if (a != req.nodes.end()) {
      samples.push_back(a->second);
    }
    if (b != req.nodes.end()) {
      samples.push_back(b->second);
    }
    if (samples.size() >= 10) {
      break;
    }
  }
  if (samples.empty()) {
    return -1;
  }

  std::set<int> preferred(preferred_faces.begin(), preferred_faces.end());
  double best_score = std::numeric_limits<double>::infinity();
  int best_idx = -1;
  int best_overlap = -1;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (edges[i].degenerated) {
      continue;
    }
    double mean = 0.0;
    int count = 0;
    for (const Vec3 &p : samples) {
      const double d = point_curve_distance(p, edges[i].edge);
      if (std::isfinite(d)) {
        mean += d;
        ++count;
      }
    }
    if (count == 0) {
      continue;
    }
    mean /= static_cast<double>(count);
    int overlap = 0;
    for (int f : edges[i].face_ids) {
      if (preferred.count(f) != 0) {
        ++overlap;
      }
    }
    if (mean < best_score ||
        (std::abs(mean - best_score) <= 1.0e-14 &&
         overlap > best_overlap)) {
      best_score = mean;
      best_idx = static_cast<int>(i);
      best_overlap = overlap;
    }
  }
  return best_idx;
}

int surface_kind(const TopoDS_Face &face) {
  try {
    BRepAdaptor_Surface s(face, Standard_False);
    return static_cast<int>(s.GetType());
  } catch (...) {
    return -1;
  }
}

bool same_support(const TopoDS_Face &a, const TopoDS_Face &b) {
  const int ka = surface_kind(a);
  const int kb = surface_kind(b);
  if (ka < 0 || ka != kb) {
    return false;
  }
  try {
    BRepAdaptor_Surface sa(a, Standard_False);
    BRepAdaptor_Surface sb(b, Standard_False);
    if (ka == 0) {
      const auto pa = sa.Plane();
      const auto pb = sb.Plane();
      const Vec3 na{pa.Axis().Direction().X(), pa.Axis().Direction().Y(),
                    pa.Axis().Direction().Z()};
      const Vec3 nb{pb.Axis().Direction().X(), pb.Axis().Direction().Y(),
                    pb.Axis().Direction().Z()};
      return norm(cross(na, nb)) < 1.0e-10;
    }
    if (ka == 1) {
      const auto ca = sa.Cylinder();
      const auto cb = sb.Cylinder();
      return std::abs(ca.Radius() - cb.Radius()) < 1.0e-8;
    }
    if (ka == 3) {
      const auto spa = sa.Sphere();
      const auto spb = sb.Sphere();
      return std::abs(spa.Radius() - spb.Radius()) < 1.0e-8;
    }
  } catch (...) {
  }
  return false;
}

std::vector<int> dedup_support_faces(const std::vector<int> &face_ids,
                                     const std::vector<OccFace> &faces) {
  std::vector<int> out;
  for (int id : face_ids) {
    if (id < 0 || id >= static_cast<int>(faces.size())) {
      continue;
    }
    bool duplicate = false;
    for (int prev : out) {
      if (same_support(faces[static_cast<std::size_t>(id)].face,
                       faces[static_cast<std::size_t>(prev)].face)) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) {
      out.push_back(id);
    }
  }
  return out;
}

bool project_surface(const TopoDS_Face &face, const Vec3 &x, Vec2 &uv,
                     Vec3 &out) {
  try {
    Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
    if (surf.IsNull()) {
      return false;
    }
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
    if (!proj.IsDone() || proj.NbPoints() <= 0) {
      return false;
    }
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

Vec3 surface_point(BRepAdaptor_Surface &surface, double u, double v) {
  gp_Pnt p = surface.Value(u, v);
  return pnt_to_vec(p);
}

std::vector<double> continuity_breaks(BRepAdaptor_Surface &surface,
                                      bool u_direction) {
  int n_intervals = 0;
  if (u_direction) {
    n_intervals = surface.NbUIntervals(GeomAbs_CN);
  } else {
    n_intervals = surface.NbVIntervals(GeomAbs_CN);
  }
  if (n_intervals <= 1) {
    return {};
  }

  TColStd_Array1OfReal intervals(1, n_intervals + 1);
  if (u_direction) {
    surface.UIntervals(intervals, GeomAbs_CN);
  } else {
    surface.VIntervals(intervals, GeomAbs_CN);
  }

  std::vector<double> values;
  values.reserve(static_cast<std::size_t>(n_intervals + 1));
  for (Standard_Integer i = intervals.Lower(); i <= intervals.Upper(); ++i) {
    values.push_back(intervals.Value(i));
  }
  if (values.size() <= 2) {
    return {};
  }

  const double first = values.front();
  const double last = values.back();
  const double tol =
      1.0e-12 * std::max({1.0, std::abs(first), std::abs(last)});

  std::vector<double> out;
  for (std::size_t i = 1; i + 1 < values.size(); ++i) {
    const double val = values[i];
    if (std::abs(val - first) <= tol || std::abs(val - last) <= tol) {
      continue;
    }
    if (!out.empty() && std::abs(val - out.back()) <= tol) {
      continue;
    }
    out.push_back(val);
  }
  return out;
}

double polyline_length(const std::vector<Vec3> &points) {
  if (points.size() < 2) {
    return 0.0;
  }
  double length = 0.0;
  for (std::size_t i = 1; i < points.size(); ++i) {
    length += norm(subtract(points[i], points[i - 1]));
  }
  return length;
}

std::vector<Vec3> deduplicate_points(const std::vector<Vec3> &points,
                                     double tol) {
  std::vector<Vec3> out;
  for (const Vec3 &p : points) {
    if (out.empty() || norm(subtract(p, out.back())) > tol) {
      out.push_back(p);
    }
  }
  return out;
}

std::vector<Vec3>
resample_polyline_by_spacing(const std::vector<Vec3> &points,
                             double target_spacing, int min_points,
                             const std::vector<Vec3> &protected_points) {
  if (points.size() <= static_cast<std::size_t>(min_points) ||
      target_spacing <= 0.0) {
    return points;
  }

  std::vector<double> seg_lengths;
  seg_lengths.reserve(points.size() - 1);
  std::vector<double> cumulative(points.size(), 0.0);
  for (std::size_t i = 1; i < points.size(); ++i) {
    const double len = norm(subtract(points[i], points[i - 1]));
    seg_lengths.push_back(len);
    cumulative[i] = cumulative[i - 1] + len;
  }

  const double total_length = cumulative.back();
  if (total_length <= 1.0e-14) {
    return points;
  }

  int n_target = static_cast<int>(std::ceil(total_length / target_spacing)) + 1;
  n_target = std::max(min_points,
                      std::min(n_target, static_cast<int>(points.size())));
  if (n_target >= static_cast<int>(points.size())) {
    return points;
  }

  std::vector<double> protected_s;
  const double projection_tol = std::max(1.0e-8, 1.0e-3 * total_length);
  for (const Vec3 &protected_point : protected_points) {
    double best_dist = std::numeric_limits<double>::infinity();
    double best_s = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t k = 0; k + 1 < points.size(); ++k) {
      const Vec3 ab = subtract(points[k + 1], points[k]);
      const double denom = dot(ab, ab);
      double alpha = 0.0;
      Vec3 q = points[k];
      if (denom > 1.0e-30) {
        alpha = dot(subtract(protected_point, points[k]), ab) / denom;
        alpha = std::max(0.0, std::min(1.0, alpha));
        q = add(points[k], scale(ab, alpha));
      }
      const double dist = norm(subtract(protected_point, q));
      if (dist < best_dist) {
        best_dist = dist;
        best_s = cumulative[k] + alpha * seg_lengths[k];
      }
    }
    if (std::isfinite(best_s) && best_dist <= projection_tol &&
        best_s > 1.0e-12 && best_s < total_length - 1.0e-12) {
      protected_s.push_back(best_s);
    }
  }

  std::vector<double> targets{0.0, total_length};
  targets.insert(targets.end(), protected_s.begin(), protected_s.end());
  const double keep_distance = 0.35 * target_spacing;
  for (int i = 1; i + 1 < n_target; ++i) {
    const double s = total_length * static_cast<double>(i) /
                     static_cast<double>(n_target - 1);
    bool near_protected = false;
    for (double ps : protected_s) {
      if (std::abs(s - ps) < keep_distance) {
        near_protected = true;
        break;
      }
    }
    if (!near_protected) {
      targets.push_back(s);
    }
  }

  std::sort(targets.begin(), targets.end());
  std::vector<double> filtered;
  const double min_sep = std::max(1.0e-12, 0.1 * target_spacing);
  for (double s : targets) {
    if (filtered.empty() || s - filtered.back() >= min_sep ||
        std::abs(s - total_length) <= 1.0e-12) {
      filtered.push_back(s);
    }
  }

  std::vector<Vec3> out;
  out.reserve(filtered.size());
  std::size_t seg_idx = 0;
  for (double s : filtered) {
    while (seg_idx + 1 < cumulative.size() - 1 &&
           cumulative[seg_idx + 1] < s) {
      ++seg_idx;
    }
    const double denom = seg_lengths[seg_idx];
    if (denom <= 1.0e-14) {
      out.push_back(points[seg_idx]);
    } else {
      const double alpha = (s - cumulative[seg_idx]) / denom;
      out.push_back(
          add(scale(points[seg_idx], 1.0 - alpha),
              scale(points[seg_idx + 1], alpha)));
    }
  }
  return out;
}

std::vector<std::vector<Vec3>>
sample_iso_segments(const TopoDS_Face &face, BRepAdaptor_Surface &surface,
                    BRepTopAdaptor_FClass2d &classifier, bool const_u,
                    double const_val, int n_samples) {
  const double umin = surface.FirstUParameter();
  const double umax = surface.LastUParameter();
  const double vmin = surface.FirstVParameter();
  const double vmax = surface.LastVParameter();
  if (!std::isfinite(umin) || !std::isfinite(umax) || !std::isfinite(vmin) ||
      !std::isfinite(vmax)) {
    return {};
  }

  const double t0 = const_u ? vmin : umin;
  const double t1 = const_u ? vmax : umax;
  if (std::abs(t1 - t0) <= 1.0e-15) {
    return {};
  }

  std::vector<std::vector<Vec3>> segments;
  std::vector<Vec3> current;
  const double point_tol = 1.0e-11;
  const int n = std::max(3, n_samples);
  for (int i = 0; i < n; ++i) {
    const double alpha = static_cast<double>(i) / static_cast<double>(n - 1);
    const double t = (1.0 - alpha) * t0 + alpha * t1;
    const double u = const_u ? const_val : t;
    const double v = const_u ? t : const_val;
    const TopAbs_State state =
        classifier.Perform(gp_Pnt2d(u, v), Standard_True);
    const bool inside = state == TopAbs_IN || state == TopAbs_ON;
    if (inside) {
      const Vec3 p = surface_point(surface, u, v);
      if (!current.empty() && norm(subtract(p, current.back())) <= point_tol) {
        continue;
      }
      current.push_back(p);
    } else {
      if (current.size() >= 2) {
        segments.push_back(current);
      }
      current.clear();
    }
  }
  if (current.size() >= 2) {
    segments.push_back(current);
  }
  (void)face;
  return segments;
}

std::vector<HelperCurve> build_face_helper_curves(const TopoDS_Face &face,
                                                  int face_tag,
                                                  double target_spacing,
                                                  int &break_count) {
  BRepAdaptor_Surface surface(face, Standard_True);
  try {
    if ((surface.IsUPeriodic() && surface.IsVPeriodic()) ||
        (surface.IsUClosed() && surface.IsVClosed())) {
      return {};
    }
  } catch (...) {
  }

  const double tol = std::max(surface.Tolerance(), 1.0e-9);
  BRepTopAdaptor_FClass2d classifier(face, tol);

  const auto u_breaks = continuity_breaks(surface, true);
  const auto v_breaks = continuity_breaks(surface, false);
  break_count += static_cast<int>(u_breaks.size() + v_breaks.size());
  if (u_breaks.empty() && v_breaks.empty()) {
    return {};
  }

  const int n_u_intervals = std::max(1, surface.NbUIntervals(GeomAbs_CN));
  const int n_v_intervals = std::max(1, surface.NbVIntervals(GeomAbs_CN));

  std::vector<HelperCurve> helpers;
  const auto add_segments = [&](const std::vector<std::vector<Vec3>> &segments,
                                const std::vector<Vec3> &protected_points,
                                int axis, double param) {
    for (const auto &segment : segments) {
      std::vector<Vec3> sampled = segment;
      if (target_spacing > 0.0) {
        sampled = resample_polyline_by_spacing(segment, target_spacing, 5,
                                               protected_points);
      }
      sampled = deduplicate_points(sampled, 1.0e-11);
      if (sampled.size() >= 2 && polyline_length(sampled) > 1.0e-12) {
        helpers.push_back({face_tag, axis, param, sampled});
      }
    }
  };

  for (double u0 : u_breaks) {
    std::vector<Vec3> protected_points;
    for (double v0 : v_breaks) {
      protected_points.push_back(surface_point(surface, u0, v0));
    }
    add_segments(sample_iso_segments(face, surface, classifier, true, u0,
                                     std::max(49, 12 * n_v_intervals + 1)),
                 protected_points, 0, u0);
  }
  for (double v0 : v_breaks) {
    std::vector<Vec3> protected_points;
    for (double u0 : u_breaks) {
      protected_points.push_back(surface_point(surface, u0, v0));
    }
    add_segments(sample_iso_segments(face, surface, classifier, false, v0,
                                     std::max(49, 12 * n_u_intervals + 1)),
                 protected_points, 1, v0);
  }
  return helpers;
}

bool surface_eval(const TopoDS_Face &face, const Vec2 &uv, SurfaceEval &e) {
  try {
    BRepAdaptor_Surface s(face, Standard_False);
    gp_Pnt p;
    gp_Vec du;
    gp_Vec dv;
    gp_Vec duu;
    gp_Vec dvv;
    gp_Vec duv;
    s.D2(uv[0], uv[1], p, du, dv, duu, dvv, duv);
    e.x = pnt_to_vec(p);
    e.xu = vec_to_vec(du);
    e.xv = vec_to_vec(dv);
    e.xuu = vec_to_vec(duu);
    e.xvv = vec_to_vec(dvv);
    e.xuv = vec_to_vec(duv);
    return finite(e.x) && finite(e.xu) && finite(e.xv);
  } catch (...) {
    return false;
  }
}

bool solve4(double a[4][4], double b[4], double x[4]) {
  double aug[4][5]{};
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      aug[i][j] = a[i][j];
    }
    aug[i][4] = b[i];
  }
  for (int c = 0; c < 4; ++c) {
    int piv = c;
    double pabs = std::abs(aug[c][c]);
    for (int r = c + 1; r < 4; ++r) {
      if (std::abs(aug[r][c]) > pabs) {
        piv = r;
        pabs = std::abs(aug[r][c]);
      }
    }
    if (pabs < 1.0e-30 || !std::isfinite(pabs)) {
      return false;
    }
    if (piv != c) {
      for (int j = c; j < 5; ++j) {
        std::swap(aug[c][j], aug[piv][j]);
      }
    }
    const double diag = aug[c][c];
    for (int j = c; j < 5; ++j) {
      aug[c][j] /= diag;
    }
    for (int r = 0; r < 4; ++r) {
      if (r == c) {
        continue;
      }
      const double f = aug[r][c];
      for (int j = c; j < 5; ++j) {
        aug[r][j] -= f * aug[c][j];
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

int normal_rank(const TopoDS_Face &f1, const TopoDS_Face &f2, const Vec3 &x) {
  Vec2 uv1{};
  Vec2 uv2{};
  Vec3 p{};
  SurfaceEval e1;
  SurfaceEval e2;
  if (!project_surface(f1, x, uv1, p) || !project_surface(f2, x, uv2, p) ||
      !surface_eval(f1, uv1, e1) || !surface_eval(f2, uv2, e2)) {
    return 0;
  }
  Vec3 n1 = cross(e1.xu, e1.xv);
  Vec3 n2 = cross(e2.xu, e2.xv);
  const double n1n = norm(n1);
  const double n2n = norm(n2);
  if (n1n < 1.0e-14 || n2n < 1.0e-14) {
    return 0;
  }
  n1 = scale(n1, 1.0 / n1n);
  n2 = scale(n2, 1.0 / n2n);
  return norm(cross(n1, n2)) > 1.0e-8 ? 2 : 1;
}

bool edge_newton(const TopoDS_Face &f1, const TopoDS_Face &f2, const Vec3 &x0,
                 Vec3 &out, double tol = 1.0e-14, int maxiter = 25) {
  Vec2 uv1{};
  Vec2 uv2{};
  Vec3 p1{};
  Vec3 p2{};
  if (!project_surface(f1, x0, uv1, p1) || !project_surface(f2, x0, uv2, p2)) {
    return false;
  }
  double y[4]{uv1[0], uv1[1], uv2[0], uv2[1]};
  SurfaceEval e1;
  SurfaceEval e2;
  for (int iter = 0; iter < maxiter; ++iter) {
    if (!surface_eval(f1, {y[0], y[1]}, e1) ||
        !surface_eval(f2, {y[2], y[3]}, e2)) {
      return false;
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
    if (!solve4(jac, rhs, dy)) {
      return false;
    }
    double step = 0.0;
    for (double v : dy) {
      step += v * v;
    }
    step = std::sqrt(step);
    if (step > 2.0) {
      for (double &v : dy) {
        v *= 2.0 / step;
      }
      step = 2.0;
    }
    for (int k = 0; k < 4; ++k) {
      y[k] += dy[k];
      if (!std::isfinite(y[k])) {
        return false;
      }
    }
    double res = 0.0;
    for (double v : r) {
      res += v * v;
    }
    if (std::sqrt(res) < tol && step < tol) {
      break;
    }
  }
  if (!surface_eval(f1, {y[0], y[1]}, e1) ||
      !surface_eval(f2, {y[2], y[3]}, e2)) {
    return false;
  }
  if (norm(subtract(e1.x, e2.x)) > 1.0e-8) {
    return false;
  }
  out = scale(add(e1.x, e2.x), 0.5);
  return finite(out);
}

Vec3 blend(const Vec3 &a, const Vec3 &b, double t) {
  return {(1.0 - t) * a[0] + t * b[0],
          (1.0 - t) * a[1] + t * b[1],
          (1.0 - t) * a[2] + t * b[2]};
}

double median_value(std::vector<double> values) {
  if (values.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::sort(values.begin(), values.end());
  const std::size_t mid = values.size() / 2;
  if (values.size() % 2 == 1) {
    return values[mid];
  }
  return 0.5 * (values[mid - 1] + values[mid]);
}

bool project_feature_node(const FeatureNode &node,
                          const std::vector<OccFace> &occ_faces,
                          const std::unordered_map<int, int> &gm_to_occ,
                          Vec3 &projected) {
  const TopoDS_Face *face = nullptr;
  for (const FeatureNodeConstraint &constraint : node.constraints) {
    const auto face_it = gm_to_occ.find(constraint.face_tag);
    if (face_it != gm_to_occ.end()) {
      face = &occ_faces[static_cast<std::size_t>(face_it->second)].face;
      break;
    }
  }
  if (face == nullptr) {
    return false;
  }

  Vec2 uv{};
  Vec3 support_point{};
  if (!project_surface(*face, node.x, uv, support_point)) {
    return false;
  }

  std::vector<double> u_values;
  std::vector<double> v_values;
  for (const FeatureNodeConstraint &constraint : node.constraints) {
    if (!std::isfinite(constraint.param)) {
      continue;
    }
    if (constraint.axis == 0) {
      u_values.push_back(constraint.param);
    } else if (constraint.axis == 1) {
      v_values.push_back(constraint.param);
    }
  }
  const double u_med = median_value(u_values);
  const double v_med = median_value(v_values);
  if (std::isfinite(u_med)) {
    uv[0] = u_med;
  }
  if (std::isfinite(v_med)) {
    uv[1] = v_med;
  }

  try {
    BRepAdaptor_Surface surface(*face, Standard_False);
    projected = surface_point(surface, uv[0], uv[1]);
    return finite(projected);
  } catch (...) {
    return false;
  }
}


} // namespace

struct Context::Impl {
  TopoDS_Shape shape;
  std::vector<OccFace> faces;
  std::vector<OccEdge> edges;
  std::vector<FaceOccurrence> occurrences;
  std::unordered_map<int, int> gm_to_occ;
  bool native_bound = false;

  explicit Impl(const std::string &step_file) : shape(read_step(step_file)),
    faces(collect_faces(shape)), edges(collect_edges(shape)) {
    TopTools_IndexedMapOfShape face_map;
    TopExp::MapShapes(shape, TopAbs_FACE, face_map);
    for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
      const int index = face_map.FindIndex(exp.Current());
      if (index < 1) throw NativeContractError("[native_identity] missing face occurrence");
      occurrences.push_back({static_cast<std::size_t>(index - 1),
                             static_cast<int>(exp.Current().Orientation())});
    }
  }
  void require_face(int tag) const {
    if (!gm_to_occ.count(tag)) throw NativeContractError(
      "[native_identity] unknown original support face tag " + std::to_string(tag));
  }
};

Context::Context(const std::string &step_file) : impl_(new Impl(step_file)) {}
Context::~Context() = default;
Context::Context(Context &&) noexcept = default;
Context &Context::operator=(Context &&) noexcept = default;
std::size_t Context::face_count() const { return impl_->faces.size(); }
const std::vector<FaceOccurrence> &Context::face_occurrences() const { return impl_->occurrences; }

std::vector<int> Context::bind_native() {
  std::string gmsh_version, build_info;
  gmsh::option::getString("General.Version", gmsh_version);
  gmsh::option::getString("General.BuildInfo", build_info);
  const std::string expected_occ = std::to_string(OCC_VERSION_MAJOR) + "." +
    std::to_string(OCC_VERSION_MINOR) + "." + std::to_string(OCC_VERSION_MAINTENANCE);
  if (gmsh_version != "4.12.2" ||
      build_info.find("OCC version: " + expected_occ + ";") == std::string::npos)
    throw NativeContractError("[native_runtime] native import requires Gmsh 4.12.2 compiled with OCCT " + expected_occ);
  for (const char *name : {"Geometry.OCCFixDegenerated", "Geometry.OCCFixSmallEdges",
                          "Geometry.OCCFixSmallFaces", "Geometry.OCCSewFaces",
                          "Geometry.OCCMakeSolids"}) {
    double value = 0.0;
    gmsh::option::getNumber(name, value);
    if (value != 0.0) throw NativeContractError(std::string("[native_runtime] native import cannot apply external file-import healing option ") + name);
  }
  double scaling = 1.0;
  gmsh::option::getNumber("Geometry.OCCScaling", scaling);
  if (scaling != 1.0) throw NativeContractError("[native_runtime] native import requires Geometry.OCCScaling=1");
  std::string target_unit;
  gmsh::option::getString("Geometry.OCCTargetUnit", target_unit);
  if (!target_unit.empty() && target_unit != "MM")
    throw NativeContractError("[native_runtime] native import preserves helper millimetre units; incompatible Geometry.OCCTargetUnit=" + target_unit);

  if (impl_->native_bound || !impl_->gm_to_occ.empty())
    throw NativeContractError("[native_identity] context has already been bound");
  gmsh::vectorpair prior;
  gmsh::model::occ::getEntities(prior);
  if (!prior.empty()) throw NativeContractError("[native_identity] native import requires an empty OCC model");
  gmsh::vectorpair imported;
  gmsh::model::occ::importShapesNativePointer(&impl_->shape, imported, true);
  gmsh::vectorpair before;
  gmsh::model::occ::getEntities(before);
  std::sort(before.begin(), before.end());
  std::vector<int> tags;
  std::set<int> unique_tags;
  for (std::size_t i = 0; i < impl_->faces.size(); ++i) {
    gmsh::vectorpair resolved;
    gmsh::model::occ::importShapesNativePointer(&impl_->faces[i].face, resolved, true);
    if (resolved.size() != 1 || resolved[0].first != 2 || resolved[0].second <= 0)
      throw NativeContractError("[native_identity] retained face native lookup did not return exactly one face; index=" + std::to_string(i));
    gmsh::vectorpair after;
    gmsh::model::occ::getEntities(after);
    std::sort(after.begin(), after.end());
    if (after != before) throw NativeContractError("[native_identity] retained face lookup changed entity inventory; index=" + std::to_string(i));
    const int tag = resolved[0].second;
    if (!unique_tags.insert(tag).second)
      throw NativeContractError("[native_identity] distinct retained faces returned the same Gmsh tag");
    impl_->gm_to_occ.emplace(tag, static_cast<int>(i));
    tags.push_back(tag);
  }
  gmsh::vectorpair model_faces;
  gmsh::model::occ::getEntities(model_faces, 2);
  if (model_faces.size() != tags.size())
    throw NativeContractError("[native_identity] native import face coverage is incomplete");
  for (const auto &entry : model_faces)
    if (!unique_tags.count(entry.second)) throw NativeContractError("[native_identity] imported face has no retained support");
  gmsh::model::occ::synchronize();
  impl_->native_bound = true;
  return tags;
}

std::vector<FaceDescriptor> Context::adapter_face_descriptors() const {
  std::vector<FaceDescriptor> out;
  for (std::size_t i = 0; i < impl_->faces.size(); ++i)
    out.push_back(face_desc_occ(impl_->faces[i].face, static_cast<int>(i)));
  return out;
}
std::string Context::adapter_brep_snapshot() const {
  std::ostringstream out;
  BRepTools::Write(impl_->shape, out, Standard_False, Standard_False,
                   TopTools_FormatVersion_CURRENT);
  return out.str();
}
void Context::bind_adapter_face_indices(const std::unordered_map<int, int> &indices) {
  if (impl_->native_bound) throw NativeContractError("[native_identity] cannot replace a native binding");
  for (const auto &entry : indices)
    if (entry.second < 0 || static_cast<std::size_t>(entry.second) >= impl_->faces.size())
      throw NativeContractError("[native_identity] invalid adapter face index");
  impl_->gm_to_occ = indices;
}

EdgeSamples Context::edge_samples(const EdgeRequest &req) const {
  const auto &occ_faces = impl_->faces;
  const auto &occ_edges = impl_->edges;
  const auto &gm_to_occ = impl_->gm_to_occ;
  for (const auto &entry : req.curve_faces)
    for (int face_tag : entry.second) impl_->require_face(face_tag);
  EdgeSamples result;

  std::unordered_map<int, std::vector<int>> curve_occ_faces;
  for (const auto &entry : req.curve_faces) {
    std::vector<int> ids;
    for (int gm_face : entry.second) {
      const auto it = gm_to_occ.find(gm_face);
      if (it != gm_to_occ.end()) {
        ids.push_back(it->second);
      }
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    curve_occ_faces[entry.first] = dedup_support_faces(ids, occ_faces);
  }

  std::unordered_map<int, int> curve_best_edge;
  for (const auto &entry : curve_occ_faces) {
    curve_best_edge[entry.first] =
        best_occ_edge_for_curve(entry.first, req, occ_edges, entry.second);
  }

  int n_newton = 0;
  int n_curve = 0;
  int n_fallback = 0;


  for (const auto &line : req.lines) {
    const auto a_it = req.nodes.find(line.a);
    const auto b_it = req.nodes.find(line.b);
    if (a_it == req.nodes.end() || b_it == req.nodes.end()) {
      throw std::runtime_error("line endpoint missing from node table");
    }
    const Vec3 p0 = a_it->second;
    const Vec3 p1 = b_it->second;
    const int edge_idx = curve_best_edge[line.curve];
    const TopoDS_Edge *edge =
        edge_idx >= 0 ? &occ_edges[static_cast<std::size_t>(edge_idx)].edge
                      : nullptr;
    const auto faces = curve_occ_faces[line.curve];
    const TopoDS_Face *f1 =
        faces.size() == 2 ? &occ_faces[static_cast<std::size_t>(faces[0])].face
                          : nullptr;
    const TopoDS_Face *f2 =
        faces.size() == 2 ? &occ_faces[static_cast<std::size_t>(faces[1])].face
                          : nullptr;

    std::vector<Vec3> samples;
    samples.reserve(req.edge_nodes.size());
    for (double a : req.edge_nodes) {
      if (std::abs(a) < 1.0e-15) {
        samples.push_back(p0);
        continue;
      }
      if (std::abs(a - 1.0) < 1.0e-15) {
        samples.push_back(p1);
        continue;
      }
      const Vec3 seed = blend(p0, p1, a);
      Vec3 on_curve = seed;
      bool projected = false;
      if (edge != nullptr && project_curve(seed, *edge, on_curve)) {
        projected = true;
        ++n_curve;
      } else {
        ++n_fallback;
      }
      Vec3 final = on_curve;
      if (f1 != nullptr && f2 != nullptr &&
          normal_rank(*f1, *f2, on_curve) >= 2) {
        Vec3 newton{};
        if (edge_newton(*f1, *f2, on_curve, newton) &&
            norm(subtract(newton, seed)) <=
                std::max(1.0e-8, 0.75 * norm(subtract(p1, p0)))) {
          final = newton;
          ++n_newton;
        }
      }
      if (!projected && !finite(final)) {
        final = seed;
      }
      samples.push_back(final);
    }
    result.lines.push_back({line.a, line.b, std::move(samples)});
  }

  result.n_newton = n_newton;
  result.n_curve = n_curve;
  result.n_fallback = n_fallback;
  return result;
}

HelperCurves Context::helper_curves(double target_spacing, const std::vector<int> &ordered_source_tags) const {
  const auto &occ_faces = impl_->faces;
  const auto &gm_to_occ = impl_->gm_to_occ;
  for (int face_tag : ordered_source_tags) impl_->require_face(face_tag);
  int break_count = 0;
  std::vector<HelperCurve> helpers;
  for (int face_tag : ordered_source_tags) {
    const auto it = gm_to_occ.find(face_tag);
    if (it == gm_to_occ.end()) {
      continue;
    }
    const auto occ_idx = static_cast<std::size_t>(it->second);
    auto face_helpers = build_face_helper_curves(
        occ_faces[occ_idx].face, face_tag, target_spacing, break_count);
    helpers.insert(helpers.end(), face_helpers.begin(), face_helpers.end());
  }

  return {std::move(helpers), break_count};
}

FeatureNodes Context::feature_nodes(const std::vector<FeatureNode> &nodes) const {
  const auto &occ_faces = impl_->faces;
  const auto &gm_to_occ = impl_->gm_to_occ;
  for (const auto &node : nodes)
    for (const auto &constraint : node.constraints) impl_->require_face(constraint.face_tag);
  FeatureNodes result;
  int n_fallback = 0;


  for (const FeatureNode &node : nodes) {
    Vec3 projected = node.x;
    bool ok = project_feature_node(node, occ_faces, gm_to_occ, projected);
    if (!ok) {
      ++n_fallback;
      projected = node.x;
    }
    const double dist = ok ? norm(subtract(projected, node.x)) : 0.0;
    result.nodes.push_back({node.tag, projected, dist});
  }

  result.n_fallback = n_fallback;
  return result;
}

FeatureSamples Context::feature_samples(const FeatureRequest &req) const {
  const auto &occ_faces = impl_->faces;
  const auto &gm_to_occ = impl_->gm_to_occ;
  for (const auto &line : req.lines) impl_->require_face(line.face_tag);
  FeatureSamples result;
  int n_fallback = 0;


  for (const FeatureLine &line : req.lines) {
    const auto a_it = req.nodes.find(line.a);
    const auto b_it = req.nodes.find(line.b);
    if (a_it == req.nodes.end() || b_it == req.nodes.end()) {
      throw std::runtime_error("feature endpoint missing from node table");
    }

    Vec3 p0 = a_it->second;
    Vec3 p1 = b_it->second;
    Vec2 uv0{};
    Vec2 uv1{};
    const auto face_it = gm_to_occ.find(line.face_tag);
    const TopoDS_Face *face = nullptr;
    if (face_it != gm_to_occ.end()) {
      face = &occ_faces[static_cast<std::size_t>(face_it->second)].face;
      Vec3 q0{};
      Vec3 q1{};
      if (project_surface(*face, p0, uv0, q0)) {
        p0 = q0;
      }
      if (project_surface(*face, p1, uv1, q1)) {
        p1 = q1;
      }
    }

    bool use_u_axis = true;
    if (line.axis == 0 && std::isfinite(line.param)) {
      use_u_axis = true;
    } else if (line.axis == 1 && std::isfinite(line.param)) {
      use_u_axis = false;
    } else {
      use_u_axis = std::abs(uv1[0] - uv0[0]) <= std::abs(uv1[1] - uv0[1]);
    }

    std::vector<Vec3> samples;
    samples.reserve(req.edge_nodes.size());
    for (double alpha : req.edge_nodes) {
      if (std::abs(alpha) < 1.0e-15) {
        samples.push_back(p0);
        continue;
      }
      if (std::abs(alpha - 1.0) < 1.0e-15) {
        samples.push_back(p1);
        continue;
      }

      Vec3 p = blend(p0, p1, alpha);
      bool ok = false;
      if (face != nullptr) {
        SurfaceEval e{};
        Vec2 uv{};
        if (use_u_axis) {
          uv = {std::isfinite(line.param)
                    ? line.param
                    : 0.5 * (uv0[0] + uv1[0]),
                (1.0 - alpha) * uv0[1] + alpha * uv1[1]};
        } else {
          uv = {(1.0 - alpha) * uv0[0] + alpha * uv1[0],
                std::isfinite(line.param)
                    ? line.param
                    : 0.5 * (uv0[1] + uv1[1])};
        }
        if (surface_eval(*face, uv, e) && finite(e.x)) {
          p = e.x;
          ok = true;
        }
      }
      if (!ok) {
        ++n_fallback;
      }
      samples.push_back(p);
    }

    result.lines.push_back({line.a, line.b, std::move(samples)});
  }

  result.n_fallback = n_fallback;
  return result;
}

ProjectedNodes Context::project_nodes(const std::vector<ProjectNode> &nodes) const {
  const auto &occ_faces = impl_->faces;
  const auto &gm_to_occ = impl_->gm_to_occ;
  for (const auto &node : nodes) impl_->require_face(node.face_tag);
  ProjectedNodes result;
  int n_fallback = 0;


  for (const ProjectNode &node : nodes) {
    Vec3 projected = node.x;
    Vec2 uv{};
    bool ok = false;
    const auto face_it = gm_to_occ.find(node.face_tag);
    if (face_it != gm_to_occ.end()) {
      const auto occ_idx = static_cast<std::size_t>(face_it->second);
      ok = project_surface(occ_faces[occ_idx].face, node.x, uv, projected);
    }
    if (!ok) {
      ++n_fallback;
      uv = {0.0, 0.0};
      projected = node.x;
    }
    const double dist = ok ? norm(subtract(projected, node.x)) : 0.0;
    result.nodes.push_back({projected, dist, uv});
  }

  result.n_fallback = n_fallback;
  return result;
}

Measurements Context::measurements(double cad_edge_tol_fraction) const {
  const TopoDS_Shape &shape = impl_->shape;

  GProp_GProps props;
  BRepGProp::SurfaceProperties(shape, props, 1.0e-9);

  const double diag = bbox_diag(shape);
  const double area = props.Mass();
  const auto &edges = impl_->edges;
  const double edge_tol = cad_edge_tol_fraction * diag;
  double min_raw_edge = std::numeric_limits<double>::infinity();
  double min_meaningful_edge = std::numeric_limits<double>::infinity();
  int n_meaningful = 0;
  for (const OccEdge &edge : edges) {
    if (edge.degenerated) {
      continue;
    }
    double len = std::numeric_limits<double>::quiet_NaN();
    try {
      len = edge_length(edge.edge);
    } catch (...) {
      continue;
    }
    if (!std::isfinite(len) || len <= 0.0) {
      continue;
    }
    min_raw_edge = std::min(min_raw_edge, len);
    if (len > edge_tol) {
      min_meaningful_edge = std::min(min_meaningful_edge, len);
      ++n_meaningful;
    }
  }
  if (!std::isfinite(min_raw_edge)) {
    min_raw_edge = std::numeric_limits<double>::quiet_NaN();
  }
  if (!std::isfinite(min_meaningful_edge)) {
    min_meaningful_edge = std::numeric_limits<double>::quiet_NaN();
  }
  if (!std::isfinite(diag) || diag <= 0.0 || !std::isfinite(area) ||
      area <= 0.0) {
    throw std::runtime_error("non-finite OCC measurements");
  }

  return {diag, area, min_raw_edge, min_meaningful_edge, edges.size(), n_meaningful};
}

} // namespace stepmesher::occ

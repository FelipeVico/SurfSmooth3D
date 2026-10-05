#include "step_mesher/occ_context.hpp"

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


namespace {
using namespace stepmesher::occ;
using FaceDesc = FaceDescriptor;
struct Node { std::size_t tag = 0; Vec3 x{}; };
struct Request : EdgeRequest { std::string step_file; std::vector<FaceDesc> gm_faces; };
struct ProjectRequest { std::string step_file; std::vector<FaceDesc> gm_faces; std::vector<ProjectNode> nodes; };
struct HelperRequest { std::string step_file; double target_spacing = 0.0; std::vector<FaceDesc> gm_faces; };
struct AdapterFeatureRequest : FeatureRequest { std::string step_file; std::vector<FaceDesc> gm_faces; };
struct FeatureNodeRequest { std::string step_file; std::vector<FaceDesc> gm_faces; std::vector<FeatureNode> nodes; };

double norm(const Vec3 &a) {
  return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

Vec3 subtract(const Vec3 &a, const Vec3 &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
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

Request read_request(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("could not open OCC edge request");
  }

  std::string marker;
  std::getline(in, marker);
  if (marker != "STEP_MESHER_OCC_EDGES_V1") {
    throw std::runtime_error("bad OCC edge request marker");
  }

  Request req;
  std::getline(in, req.step_file);

  std::size_t n = 0;
  in >> n;
  req.edge_nodes.resize(n);
  for (double &a : req.edge_nodes) {
    in >> a;
  }

  in >> n;
  req.gm_faces.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    int tag = 0;
    std::array<double, 6> bbox{};
    double area = 0.0;
    in >> tag >> bbox[0] >> bbox[1] >> bbox[2] >> bbox[3] >> bbox[4] >>
        bbox[5] >> area;
    req.gm_faces.push_back(make_desc(bbox, area, tag));
  }

  in >> n;
  for (std::size_t i = 0; i < n; ++i) {
    Node node;
    in >> node.tag >> node.x[0] >> node.x[1] >> node.x[2];
    req.nodes[node.tag] = node.x;
  }

  in >> n;
  for (std::size_t i = 0; i < n; ++i) {
    int curve = 0;
    std::size_t k = 0;
    in >> curve >> k;
    auto &faces = req.curve_faces[curve];
    faces.resize(k);
    for (int &face : faces) {
      in >> face;
    }
  }

  in >> n;
  req.lines.resize(n);
  for (Line &line : req.lines) {
    in >> line.a >> line.b >> line.curve;
  }

  if (!in) {
    throw std::runtime_error("truncated OCC edge request");
  }
  return req;
}

ProjectRequest read_project_request(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("could not open OCC projection request");
  }

  std::string marker;
  std::getline(in, marker);
  if (marker != "STEP_MESHER_OCC_PROJECT_V1") {
    throw std::runtime_error("bad OCC projection request marker");
  }

  ProjectRequest req;
  std::getline(in, req.step_file);

  std::size_t n = 0;
  in >> n;
  req.gm_faces.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    int tag = 0;
    std::array<double, 6> bbox{};
    double area = 0.0;
    in >> tag >> bbox[0] >> bbox[1] >> bbox[2] >> bbox[3] >> bbox[4] >>
        bbox[5] >> area;
    req.gm_faces.push_back(make_desc(bbox, area, tag));
  }

  in >> n;
  req.nodes.resize(n);
  for (ProjectNode &node : req.nodes) {
    in >> node.face_tag >> node.x[0] >> node.x[1] >> node.x[2];
  }

  if (!in) {
    throw std::runtime_error("truncated OCC projection request");
  }
  return req;
}

std::vector<FaceDesc> read_face_descs(std::istream &in) {
  std::size_t n = 0;
  in >> n;
  std::vector<FaceDesc> faces;
  faces.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    int tag = 0;
    std::array<double, 6> bbox{};
    double area = 0.0;
    in >> tag >> bbox[0] >> bbox[1] >> bbox[2] >> bbox[3] >> bbox[4] >>
        bbox[5] >> area;
    faces.push_back(make_desc(bbox, area, tag));
  }
  return faces;
}

HelperRequest read_helper_request(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("could not open OCC helper-curve request");
  }

  std::string marker;
  std::getline(in, marker);
  if (marker != "STEP_MESHER_OCC_HELPER_CURVES_V1") {
    throw std::runtime_error("bad OCC helper-curve request marker");
  }

  HelperRequest req;
  std::getline(in, req.step_file);
  in >> req.target_spacing;
  req.gm_faces = read_face_descs(in);

  if (!in) {
    throw std::runtime_error("truncated OCC helper-curve request");
  }
  return req;
}

AdapterFeatureRequest read_feature_request(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("could not open OCC feature-edge request");
  }

  std::string marker;
  std::getline(in, marker);
  if (marker != "STEP_MESHER_OCC_FEATURE_EDGES_V1") {
    throw std::runtime_error("bad OCC feature-edge request marker");
  }

  AdapterFeatureRequest req;
  std::getline(in, req.step_file);

  std::size_t n = 0;
  in >> n;
  req.edge_nodes.resize(n);
  for (double &a : req.edge_nodes) {
    in >> a;
  }

  req.gm_faces = read_face_descs(in);

  in >> n;
  for (std::size_t i = 0; i < n; ++i) {
    Node node;
    in >> node.tag >> node.x[0] >> node.x[1] >> node.x[2];
    req.nodes[node.tag] = node.x;
  }

  in >> n;
  req.lines.resize(n);
  for (FeatureLine &line : req.lines) {
    in >> line.a >> line.b >> line.face_tag >> line.axis >> line.param;
  }

  if (!in) {
    throw std::runtime_error("truncated OCC feature-edge request");
  }
  return req;
}

FeatureNodeRequest read_feature_node_request(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("could not open OCC feature-node request");
  }

  std::string marker;
  std::getline(in, marker);
  if (marker != "STEP_MESHER_OCC_FEATURE_NODES_V1") {
    throw std::runtime_error("bad OCC feature-node request marker");
  }

  FeatureNodeRequest req;
  std::getline(in, req.step_file);
  req.gm_faces = read_face_descs(in);

  std::size_t n = 0;
  in >> n;
  req.nodes.resize(n);
  for (FeatureNode &node : req.nodes) {
    std::size_t k = 0;
    in >> node.tag >> node.x[0] >> node.x[1] >> node.x[2] >> k;
    node.constraints.resize(k);
    for (FeatureNodeConstraint &constraint : node.constraints) {
      in >> constraint.face_tag >> constraint.axis >> constraint.param;
    }
  }

  if (!in) {
    throw std::runtime_error("truncated OCC feature-node request");
  }
  return req;
}

double desc_cost(const FaceDesc &g, const FaceDesc &o, double model_diag,
                 double area_scale) {
  double bbox_term = 0.0;
  for (int i = 0; i < 6; ++i) {
    const double d = g.bbox[i] - o.bbox[i];
    bbox_term += d * d;
  }
  bbox_term = std::sqrt(bbox_term) / std::max(model_diag, 1.0e-15);
  const double center_term =
      norm(subtract(g.center, o.center)) / std::max(model_diag, 1.0e-15);
  const double size_term =
      norm(subtract(g.size, o.size)) / std::max(model_diag, 1.0e-15);
  const double diag_term =
      std::abs(g.diag - o.diag) / std::max(model_diag, 1.0e-15);
  const double area_term =
      std::abs(g.area - o.area) / std::max(area_scale, 1.0e-15);
  return 5.0 * bbox_term + 4.0 * center_term + 4.0 * size_term +
         2.0 * diag_term + 2.0 * area_term;
}

std::unordered_map<int, int> map_faces(const std::vector<FaceDesc> &gm_faces,
                                       const std::vector<FaceDesc> &occ_faces,
                                       double model_diag) {
  double area_scale = 1.0;
  for (const auto &f : gm_faces) {
    area_scale = std::max(area_scale, std::abs(f.area));
  }
  for (const auto &f : occ_faces) {
    area_scale = std::max(area_scale, std::abs(f.area));
  }

  std::unordered_map<int, int> out;
  std::set<int> used;
  for (const auto &gm_face : gm_faces) {
    double best = std::numeric_limits<double>::infinity();
    int best_idx = -1;
    for (std::size_t j = 0; j < occ_faces.size(); ++j) {
      if (used.count(static_cast<int>(j)) != 0) {
        continue;
      }
      const double c =
          desc_cost(gm_face, occ_faces[j], model_diag, area_scale);
      if (c < best) {
        best = c;
        best_idx = static_cast<int>(j);
      }
    }
    if (best_idx < 0) {
      for (std::size_t j = 0; j < occ_faces.size(); ++j) {
        const double c =
            desc_cost(gm_face, occ_faces[j], model_diag, area_scale);
        if (c < best) {
          best = c;
          best_idx = static_cast<int>(j);
        }
      }
    }
    if (best_idx >= 0) {
      out[gm_face.tag] = best_idx;
      used.insert(best_idx);
    }
  }
  return out;
}

Context context_for_adapter(const std::string &step_file, const std::vector<FaceDesc> &descriptors) {
  Context context(step_file);
  const auto measurements = context.measurements(1e-6);
  context.bind_adapter_face_indices(map_faces(descriptors, context.adapter_face_descriptors(), measurements.diag));
  return context;
}
void write_lines(std::ostream &out, const std::vector<LineSamples> &lines) {
  for (const auto &line : lines) {
    out << line.a << " " << line.b << " " << line.points.size();
    for (const auto &p : line.points) out << " " << p[0] << " " << p[1] << " " << p[2];
    out << "\n";
  }
}
std::ofstream output_file(const std::string &path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("could not open OCC output");
  out << std::setprecision(17);
  return out;
}
void write_edge_samples(const Request &req, const std::string &path) {
  auto context = context_for_adapter(req.step_file, req.gm_faces);
  const auto result = context.edge_samples(req);
  auto out = output_file(path);
  out << "OK " << result.lines.size() << " " << result.n_newton << " " << result.n_curve << " " << result.n_fallback << "\n";
  write_lines(out, result.lines);
}
void write_helper_curves(const HelperRequest &req, const std::string &path) {
  auto context = context_for_adapter(req.step_file, req.gm_faces);
  std::vector<int> tags;
  for (const auto &face : req.gm_faces) tags.push_back(face.tag);
  const auto result = context.helper_curves(req.target_spacing, tags);
  auto out = output_file(path);
  out << "OK " << result.curves.size() << " " << result.break_count << "\n";
  for (const auto &helper : result.curves) {
    out << helper.face_tag << " " << helper.axis << " " << helper.param << " " << helper.points.size();
    for (const auto &p : helper.points) out << " " << p[0] << " " << p[1] << " " << p[2];
    out << "\n";
  }
}
void write_feature_nodes(const FeatureNodeRequest &req, const std::string &path) {
  auto context = context_for_adapter(req.step_file, req.gm_faces);
  const auto result = context.feature_nodes(req.nodes);
  auto out = output_file(path);
  out << "OK " << result.nodes.size() << " " << result.n_fallback << "\n";
  for (const auto &node : result.nodes) out << node.tag << " " << node.x[0] << " " << node.x[1] << " " << node.x[2] << " " << node.distance << "\n";
}
void write_feature_samples(const AdapterFeatureRequest &req, const std::string &path) {
  auto context = context_for_adapter(req.step_file, req.gm_faces);
  const auto result = context.feature_samples(req);
  auto out = output_file(path);
  out << "OK " << result.lines.size() << " " << result.n_fallback << "\n";
  write_lines(out, result.lines);
}
void write_projected_nodes(const ProjectRequest &req, const std::string &path) {
  auto context = context_for_adapter(req.step_file, req.gm_faces);
  const auto result = context.project_nodes(req.nodes);
  auto out = output_file(path);
  out << "OK " << result.nodes.size() << " " << result.n_fallback << "\n";
  for (const auto &node : result.nodes) out << node.x[0] << " " << node.x[1] << " " << node.x[2] << " " << node.distance << " " << node.uv[0] << " " << node.uv[1] << "\n";
}
void write_info(const std::string &step_file, double fraction) {
  const auto result = Context(step_file).measurements(fraction);
  std::cout << std::setprecision(17) << "OK " << result.diag << " " << result.area << " " << result.min_raw_edge << " " << result.min_meaningful_edge << " " << result.n_edges << " " << result.n_meaningful << "\n";
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc == 2) {
      write_info(argv[1], 1.0e-6);
      return 0;
    }
    if (argc == 3) {
      write_info(argv[1], std::stod(argv[2]));
      return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "edge-samples") {
      write_edge_samples(read_request(argv[2]), argv[3]);
      return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "helper-curves") {
      write_helper_curves(read_helper_request(argv[2]), argv[3]);
      return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "feature-samples") {
      write_feature_samples(read_feature_request(argv[2]), argv[3]);
      return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "feature-nodes") {
      write_feature_nodes(read_feature_node_request(argv[2]), argv[3]);
      return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "project-nodes") {
      write_projected_nodes(read_project_request(argv[2]), argv[3]);
      return 0;
    }
    std::cerr << "Usage:\n"
              << "  step_mesher_occ_info file.step\n"
              << "  step_mesher_occ_info file.step cad_edge_tol_fraction\n"
              << "  step_mesher_occ_info edge-samples input.txt output.txt\n"
              << "  step_mesher_occ_info helper-curves input.txt output.txt\n"
              << "  step_mesher_occ_info feature-samples input.txt output.txt\n"
              << "  step_mesher_occ_info feature-nodes input.txt output.txt\n"
              << "  step_mesher_occ_info project-nodes input.txt output.txt\n";
    return 2;
  } catch (const std::exception &e) {
    std::cerr << "step_mesher_occ_info: " << e.what() << "\n";
    return 1;
  }
}

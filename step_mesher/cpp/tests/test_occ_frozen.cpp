// Test-only adapter reuse supplies the unchanged legacy input grammar. This is
// not linked into either production entry point.
#define main standalone_occ_adapter_main
#include "../src/occ_info.cpp"
#undef main

int main(int argc, char **argv) {
  try {
    if (argc != 5) throw std::runtime_error("usage: test_occ_frozen operation request frozen-map output; operation freeze:<operation> records a comparison-only face association");
    std::string mode = argv[1];
    bool freeze = mode.rfind("freeze:", 0) == 0;
    if (freeze) mode.erase(0, 7);
    std::string step_file;
    std::vector<FaceDesc> descriptors;
    Request edge;
    ProjectRequest projected;
    HelperRequest helper;
    AdapterFeatureRequest feature;
    FeatureNodeRequest nodes;
    if (mode == "edge-samples") { edge = read_request(argv[2]); step_file = edge.step_file; descriptors = edge.gm_faces; }
    else if (mode == "project-nodes") { projected = read_project_request(argv[2]); step_file = projected.step_file; descriptors = projected.gm_faces; }
    else if (mode == "helper-curves") { helper = read_helper_request(argv[2]); step_file = helper.step_file; descriptors = helper.gm_faces; }
    else if (mode == "feature-samples") { feature = read_feature_request(argv[2]); step_file = feature.step_file; descriptors = feature.gm_faces; }
    else if (mode == "feature-nodes") { nodes = read_feature_node_request(argv[2]); step_file = nodes.step_file; descriptors = nodes.gm_faces; }
    else throw std::runtime_error("unsupported frozen operation");
    Context context(step_file);
    std::unordered_map<int, int> indices;
    if (freeze) {
      const auto measured = context.measurements(1e-6);
      indices = map_faces(descriptors, context.adapter_face_descriptors(), measured.diag);
      std::ofstream mapping(argv[3]);
      if (!mapping) throw std::runtime_error("could not create frozen association");
      std::map<int, int> ordered(indices.begin(), indices.end());
      for (const auto &entry : ordered) mapping << entry.first << " " << entry.second << "\n";
      mapping.close();
      if (!mapping) throw std::runtime_error("could not close frozen association");
    } else {
      std::ifstream mapping(argv[3]);
      if (!mapping) throw std::runtime_error("could not open frozen association");
      int tag = 0, index = 0;
      while (mapping >> tag >> index) {
        if (!indices.emplace(tag, index).second) throw std::runtime_error("duplicate frozen tag");
      }
      if (!mapping.eof()) throw std::runtime_error("malformed frozen association");
    }
    context.bind_adapter_face_indices(indices);
    auto out = output_file(argv[4]);
    if (mode == "edge-samples") {
      const auto result = context.edge_samples(edge);
      out << "OK " << result.lines.size() << " " << result.n_newton << " " << result.n_curve << " " << result.n_fallback << "\n";
      write_lines(out, result.lines);
    } else if (mode == "feature-samples") {
      const auto result = context.feature_samples(feature);
      out << "OK " << result.lines.size() << " " << result.n_fallback << "\n";
      write_lines(out, result.lines);
    } else if (mode == "helper-curves") {
      std::vector<int> tags;
      for (const auto &face : helper.gm_faces) tags.push_back(face.tag);
      const auto result = context.helper_curves(helper.target_spacing, tags);
      out << "OK " << result.curves.size() << " " << result.break_count << "\n";
      for (const auto &curve : result.curves) {
        out << curve.face_tag << " " << curve.axis << " " << curve.param << " " << curve.points.size();
        for (const auto &point : curve.points) out << " " << point[0] << " " << point[1] << " " << point[2];
        out << "\n";
      }
    } else if (mode == "feature-nodes") {
      const auto result = context.feature_nodes(nodes.nodes);
      out << "OK " << result.nodes.size() << " " << result.n_fallback << "\n";
      for (const auto &node : result.nodes) out << node.tag << " " << node.x[0] << " " << node.x[1] << " " << node.x[2] << " " << node.distance << "\n";
    } else {
      const auto result = context.project_nodes(projected.nodes);
      out << "OK " << result.nodes.size() << " " << result.n_fallback << "\n";
      for (const auto &node : result.nodes) out << node.x[0] << " " << node.x[1] << " " << node.x[2] << " " << node.distance << " " << node.uv[0] << " " << node.uv[1] << "\n";
    }
    out.close();
    if (!out) throw std::runtime_error("could not close frozen output");
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "test_occ_frozen: " << error.what() << "\n";
    return 1;
  }
}

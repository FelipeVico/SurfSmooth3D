#define main standalone_occ_adapter_main
#include "../src/occ_info.cpp"
#undef main
#include <gmsh.h>

namespace {
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

void checked_write(const std::string &path, const std::string &text) {
  std::ofstream out(path); out << text; out.close();
  if (!out) throw std::runtime_error("could not write retention evidence");
}
}
int main(int argc,char **argv) {
  bool initialized=false;
  try {
    if(argc!=5)throw std::runtime_error("usage: test_occ_retention operation request helper-request output-prefix");
    const std::string mode=argv[1];
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

    gmsh::initialize(0,nullptr,false,false);initialized=true;
    gmsh::option::setNumber("General.Terminal",0);
    gmsh::model::add("retention-fragment-test");
    Context context(step_file);
    const auto tags=context.bind_native();
    const auto measures=context.measurements(1e-6);
    const auto old_map=map_faces(descriptors,context.adapter_face_descriptors(),measures.diag);
    auto remap=[&](int tag) { return tags.at(static_cast<std::size_t>(old_map.at(tag))); };
    for(auto &entry:edge.curve_faces)for(int &face:entry.second)face=remap(face);
    for(auto &node:projected.nodes)node.face_tag=remap(node.face_tag);
    for(auto &line:feature.lines)line.face_tag=remap(line.face_tag);
    for(auto &node:nodes.nodes)for(auto &constraint:node.constraints)constraint.face_tag=remap(constraint.face_tag);
    for(auto &face:helper.gm_faces)face.tag=remap(face.tag);
    auto evaluate=[&]() {
    std::ostringstream out;
    out << std::setprecision(17);
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

      return out.str();
    };
    const auto before_result=evaluate();
    const auto before_brep=context.adapter_brep_snapshot();
    gmsh::vectorpair volumes;gmsh::model::getEntities(volumes,3);
    if(!volumes.empty())gmsh::model::occ::remove(volumes,false);
    gmsh::model::occ::synchronize();
    const auto helper_request=read_helper_request(argv[3]);
    std::vector<int> ordered_tags;
    for(const auto &face:helper_request.gm_faces)ordered_tags.push_back(remap(face.tag));
    const auto curves=context.helper_curves(helper_request.target_spacing,ordered_tags);
    gmsh::vectorpair tool_curves,object_faces;
    for(const auto &curve:curves.curves) {
      const int tag=add_helper_curve_to_gmsh(curve.points);
      if(tag>0)tool_curves.emplace_back(1,tag);
    }
    for(int tag:tags)object_faces.emplace_back(2,tag);
    if(!tool_curves.empty()) {
      gmsh::model::occ::synchronize();
      gmsh::vectorpair outputs;std::vector<gmsh::vectorpair> ancestry;
      gmsh::model::occ::fragment(object_faces,tool_curves,outputs,ancestry,-1,true,true);
      gmsh::model::occ::synchronize();
    }
    const auto after_result=evaluate();
    const auto after_brep=context.adapter_brep_snapshot();
    const auto after_measures=context.measurements(1e-6);
    const bool measures_equal=measures.diag==after_measures.diag && measures.area==after_measures.area &&
      measures.min_raw_edge==after_measures.min_raw_edge && measures.min_meaningful_edge==after_measures.min_meaningful_edge &&
      measures.n_edges==after_measures.n_edges && measures.n_meaningful==after_measures.n_meaningful;
    const std::string prefix=argv[4];
    checked_write(prefix+".before.txt",before_result);checked_write(prefix+".after.txt",after_result);
    checked_write(prefix+".before.brep",before_brep);checked_write(prefix+".after.brep",after_brep);
    std::cout<<"{\"operations_byte_identical\":"<<(before_result==after_result?"true":"false")<<",\"brep_byte_identical\":"<<(before_brep==after_brep?"true":"false")<<",\"measurements_equal\":"<<(measures_equal?"true":"false")<<",\"helpers\":"<<tool_curves.size()<<"}\n";
    gmsh::finalize();return before_result==after_result?0:1;
  }catch(const std::exception &error){if(initialized)gmsh::finalize();std::cerr<<error.what()<<"\n";return 1;}
  catch(...){if(initialized)gmsh::finalize();std::cerr<<"unknown retention failure\n";return 1;}
}

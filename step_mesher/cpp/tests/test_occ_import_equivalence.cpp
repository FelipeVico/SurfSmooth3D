// Descriptor association is deliberately confined to this comparison test for
// the old independent Gmsh file import. Production uses native topology only.
#define main standalone_occ_adapter_main
#include "../src/occ_info.cpp"
#undef main
#include <gmsh.h>

namespace {
struct SurfaceSnapshot {
  std::vector<double> lower, upper, xyz, first, second;
  std::string type;
};
struct ImportSnapshot {
  std::vector<SurfaceSnapshot> surfaces;
  std::map<std::vector<int>, std::vector<double>> edge_lengths_by_adjacency;
  std::vector<std::pair<int,int>> original_association;
};
std::vector<FaceDescriptor> gmsh_descriptors() {
  gmsh::vectorpair entities;
  gmsh::model::getEntities(entities, 2);
  std::vector<FaceDescriptor> result;
  for (const auto &entry : entities) {
    std::array<double,6> bounds{};
    gmsh::model::occ::getBoundingBox(2, entry.second, bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]);
    double area = 0;
    gmsh::model::occ::getMass(2,entry.second,area);
    result.push_back(make_desc(bounds,area,entry.second));
  }
  return result;
}
ImportSnapshot snapshot(const std::unordered_map<int,int> &mapping, std::size_t count) {
  ImportSnapshot result;
  result.surfaces.resize(count);
  std::set<int> indices;
  for (const auto &entry : mapping) {
    if (entry.second < 0 || static_cast<std::size_t>(entry.second) >= count || !indices.insert(entry.second).second)
      throw std::runtime_error("comparison association is not one-to-one");
    result.original_association.push_back(entry);
    auto &surface = result.surfaces[static_cast<std::size_t>(entry.second)];
    gmsh::model::getType(2,entry.first,surface.type);
    gmsh::model::getParametrizationBounds(2,entry.first,surface.lower,surface.upper);
    if (surface.lower.size()!=2 || surface.upper.size()!=2) throw std::runtime_error("invalid face parameter bounds");
    std::vector<double> uv;
    for(int j=0;j<9;++j)for(int i=0;i<9;++i) {
      uv.push_back(surface.lower[0]+(surface.upper[0]-surface.lower[0])*(double(i)/8.0));
      uv.push_back(surface.lower[1]+(surface.upper[1]-surface.lower[1])*(double(j)/8.0));
    }
    gmsh::model::getValue(2,entry.first,uv,surface.xyz);
    gmsh::model::getDerivative(2,entry.first,uv,surface.first);
    gmsh::model::getSecondDerivative(2,entry.first,uv,surface.second);
  }
  if(indices.size()!=count)throw std::runtime_error("comparison association incomplete");
  std::sort(result.original_association.begin(),result.original_association.end());
  gmsh::vectorpair edges;
  gmsh::model::getEntities(edges,1);
  for(const auto &edge:edges) {
    std::vector<int> upwards,downwards,faces;
    gmsh::model::getAdjacencies(1,edge.second,upwards,downwards);
    for(int face:upwards)faces.push_back(mapping.at(face));
    std::sort(faces.begin(),faces.end());
    double length=0;
    gmsh::model::occ::getMass(1,edge.second,length);
    result.edge_lengths_by_adjacency[faces].push_back(length);
  }
  for(auto &entry:result.edge_lengths_by_adjacency)std::sort(entry.second.begin(),entry.second.end());
  return result;
}
void compare_vectors(const std::vector<double> &a,const std::vector<double> &b,double &worst,double &scale) {
  if(a.size()!=b.size())throw std::runtime_error("different evaluation sizes");
  for(std::size_t i=0;i<a.size();++i){
    if(!std::isfinite(a[i])||!std::isfinite(b[i]))throw std::runtime_error("nonfinite support evaluation");
    worst=std::max(worst,std::abs(a[i]-b[i]));scale=std::max({scale,std::abs(a[i]),std::abs(b[i])});
  }
}
}
int main(int argc,char **argv) {
  if(argc!=2)return 2;
  bool initialized=false;
  try {
    gmsh::initialize(0,nullptr,false,false);initialized=true;
    gmsh::option::setNumber("General.Terminal",0);
    gmsh::model::add("native-retained");
    Context context(argv[1]);
    const auto measure=context.measurements(1e-6);
    const auto retained=context.adapter_face_descriptors();
    const auto tags=context.bind_native();
    std::unordered_map<int,int> native_map;
    for(std::size_t i=0;i<tags.size();++i)native_map.emplace(tags[i],static_cast<int>(i));
    const auto native=snapshot(native_map,tags.size());
    gmsh::clear();gmsh::model::add("old-file-reader");
    gmsh::vectorpair imported;
    gmsh::model::occ::importShapes(argv[1],imported,true);
    gmsh::model::occ::synchronize();
    const auto old_desc=gmsh_descriptors();
    if(old_desc.size()!=tags.size())throw std::runtime_error("file import and retained import have different face counts");
    const auto old_map=map_faces(old_desc,retained,measure.diag);
    const auto old=snapshot(old_map,tags.size());
    double position_error=0,position_scale=measure.diag,derivative_error=0,derivative_scale=0,bounds_error=0,bounds_scale=0,edge_error=0,edge_scale=measure.diag;
    for(std::size_t i=0;i<tags.size();++i) {
      const auto &a=native.surfaces[i],&b=old.surfaces[i];
      if(a.type!=b.type)throw std::runtime_error("different support types at source face "+std::to_string(i));
      compare_vectors(a.lower,b.lower,bounds_error,bounds_scale);
      compare_vectors(a.upper,b.upper,bounds_error,bounds_scale);
      compare_vectors(a.xyz,b.xyz,position_error,position_scale);
      compare_vectors(a.first,b.first,derivative_error,derivative_scale);
      compare_vectors(a.second,b.second,derivative_error,derivative_scale);
    }
    if(native.edge_lengths_by_adjacency.size()!=old.edge_lengths_by_adjacency.size())throw std::runtime_error("different source face adjacency");
    for(const auto &entry:native.edge_lengths_by_adjacency) {
      const auto found=old.edge_lengths_by_adjacency.find(entry.first);
      if(found==old.edge_lengths_by_adjacency.end())throw std::runtime_error("missing source face adjacency group");
      compare_vectors(entry.second,found->second,edge_error,edge_scale);
    }
    const double tolerance=128*std::numeric_limits<double>::epsilon();
    const bool passed=position_error<=tolerance*position_scale&&derivative_error<=tolerance*derivative_scale&&bounds_error<=tolerance*bounds_scale&&edge_error<=tolerance*edge_scale;
    std::cout<<std::setprecision(17)<<"{\"passed\":"<<(passed?"true":"false")<<",\"faces\":"<<tags.size()<<",\"samples_per_support\":81,\"max_position_difference\":"<<position_error<<",\"max_derivative_difference\":"<<derivative_error<<",\"max_bound_difference\":"<<bounds_error<<",\"max_edge_length_difference\":"<<edge_error<<",\"position_limit\":"<<tolerance*position_scale<<",\"same_numeric_tag_association\":"<<(native.original_association==old.original_association?"true":"false")<<",\"native_association\":[";
    for(std::size_t i=0;i<native.original_association.size();++i){if(i)std::cout<<",";const auto &p=native.original_association[i];std::cout<<"["<<p.first<<","<<p.second<<"]";}
    std::cout<<"],\"old_association\":[";
    for(std::size_t i=0;i<old.original_association.size();++i){if(i)std::cout<<",";const auto &p=old.original_association[i];std::cout<<"["<<p.first<<","<<p.second<<"]";}
    std::cout<<"]}\n";
    gmsh::finalize();return passed?0:1;
  }catch(const std::exception &error){if(initialized)gmsh::finalize();std::cerr<<error.what()<<"\n";return 1;}
  catch(...){if(initialized)gmsh::finalize();std::cerr<<"unknown import probe failure\n";return 1;}
}

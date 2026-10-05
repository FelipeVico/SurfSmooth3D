#include "occ_context_probe.hpp"
#include "step_mesher/occ_context.hpp"
#include <gmsh.h>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Plane.hxx>
#include <Precision.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRep_Builder.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Controller.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Trsf.hxx>
#include <gp_Pln.hxx>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
void require(bool ok, const std::string &message) {
  if (!ok) throw std::runtime_error("native OCC probe: " + message);
}
gmsh::vectorpair inventory() {
  gmsh::vectorpair out;
  gmsh::model::occ::getEntities(out);
  std::sort(out.begin(), out.end());
  return out;
}
void retain_entity_check(const TopoDS_Shape &shape, int dimension, TopAbs_ShapeEnum kind,
                         const gmsh::vectorpair &before) {
  TopTools_IndexedMapOfShape entities;
  TopExp::MapShapes(shape, kind, entities);
  std::set<int> tags;
  for (int i = 1; i <= entities.Extent(); ++i) {
    const TopoDS_Shape &entity = entities.FindKey(i);
    gmsh::vectorpair resolved;
    gmsh::model::occ::importShapesNativePointer(&entity, resolved, true);
    require(resolved.size() == 1 && resolved[0].first == dimension,
            "native lookup must return one existing entity");
    require(tags.insert(resolved[0].second).second, "distinct topology must retain distinct tags");
    require(inventory() == before, "lookup must not increase entity counts");
    const TopoDS_Shape reverse = entity.Reversed();
    gmsh::vectorpair reversed;
    gmsh::model::occ::importShapesNativePointer(&reverse, reversed, true);
    require(reversed == resolved && inventory() == before,
            "orientation must not change canonical identity");
  }
}
void native_fixture(bool coincident) {
  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(1.0, 2.0, 3.0).Shape();
  if (!coincident) {
    builder.Add(compound, box);
    gp_Trsf transform;
    transform.SetTranslation(gp_Vec(4.0, 0.0, 0.0));
    builder.Add(compound, box.Moved(TopLoc_Location(transform)));
    TopExp_Explorer face(box, TopAbs_FACE);
    builder.Add(compound, face.Current().Reversed());
  } else {
    TopExp_Explorer face(box, TopAbs_FACE);
    builder.Add(compound, face.Current());
    builder.Add(compound, BRepBuilderAPI_Copy(face.Current()).Shape());
  }
  gmsh::vectorpair imported;
  gmsh::model::occ::importShapesNativePointer(&compound, imported, true);
  const auto before = inventory();
  retain_entity_check(compound, 2, TopAbs_FACE, before);
  retain_entity_check(compound, 1, TopAbs_EDGE, before);
  retain_entity_check(compound, 0, TopAbs_VERTEX, before);
  if (!coincident) {
    gmsh::model::occ::synchronize();
    gmsh::option::setNumber("Mesh.MeshSizeMin", 0.7);
    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.7);
    gmsh::model::mesh::generate(2);
    std::vector<std::size_t> tags;
    std::vector<double> xyz, uv;
    gmsh::model::mesh::getNodes(tags, xyz, uv);
    require(!tags.empty(), "retained located shape must mesh");
  }
}
void shared_support_trim_fixture() {
  Handle(Geom_Surface) support = new Geom_Plane(gp_Pln(gp_Pnt(0,0,0),gp_Dir(0,0,1)));
  const TopoDS_Face square = BRepBuilderAPI_MakeFace(support,-1,1,-1,1,Precision::Confusion()).Face();
  BRepBuilderAPI_MakePolygon polygon;
  polygon.Add(gp_Pnt(-1,-1,0)); polygon.Add(gp_Pnt(1,-1,0)); polygon.Add(gp_Pnt(-1,1,0)); polygon.Close();
  const TopoDS_Face triangle = BRepBuilderAPI_MakeFace(support,polygon.Wire(),Standard_True).Face();
  require(BRep_Tool::Surface(square) == support && BRep_Tool::Surface(triangle) == support,
          "trim fixture must share the exact original support handle");
  require(!square.IsSame(triangle), "different trims must be distinct topology");
  BRep_Builder builder; TopoDS_Compound compound; builder.MakeCompound(compound);
  builder.Add(compound,square); builder.Add(compound,triangle);
  gmsh::vectorpair imported;
  gmsh::model::occ::importShapesNativePointer(&compound,imported,true);
  require(imported.size()==2 && imported[0].second!=imported[1].second,
          "shared support and equal bounding boxes must not merge distinct trims");
  const auto before=inventory();
  retain_entity_check(compound,2,TopAbs_FACE,before);
  retain_entity_check(compound,1,TopAbs_EDGE,before);
  retain_entity_check(compound,0,TopAbs_VERTEX,before);
}
void periodic_fixture() {
  const TopoDS_Shape torus=BRepPrimAPI_MakeTorus(3.0,1.0).Shape();
  gmsh::vectorpair imported;
  gmsh::model::occ::importShapesNativePointer(&torus,imported,true);
  const auto before=inventory();
  retain_entity_check(torus,2,TopAbs_FACE,before);
  retain_entity_check(torus,1,TopAbs_EDGE,before);
  retain_entity_check(torus,0,TopAbs_VERTEX,before);
  gmsh::model::occ::synchronize();
  gmsh::option::setNumber("Mesh.MeshSizeMin",0.7);
  gmsh::option::setNumber("Mesh.MeshSizeMax",0.7);
  gmsh::model::mesh::generate(2);
}
void fragment_fixture() {
  const int face=gmsh::model::occ::addRectangle(0,0,0,2,2);
  const int first=gmsh::model::occ::addPoint(1,0,0);
  const int last=gmsh::model::occ::addPoint(1,2,0);
  const int line=gmsh::model::occ::addLine(first,last);
  gmsh::model::occ::synchronize();
  gmsh::vectorpair result;std::vector<gmsh::vectorpair> ancestry;
  gmsh::model::occ::fragment({{2,face}},{{1,line}},result,ancestry,-1,true,true);
  gmsh::model::occ::synchronize();
  require(ancestry.size()==2,"fragment must return object and tool ancestry");
  std::set<int> child_faces;
  for(const auto &entry:ancestry[0])if(entry.first==2)child_faces.insert(entry.second);
  require(child_faces.size()==2,"partition line must produce two face descendants");
  gmsh::vectorpair live;gmsh::model::getEntities(live,2);
  require(live.size()==child_faces.size(),"every live face must be covered by parent ancestry");
  for(const auto &entry:live) {
    require(child_faces.count(entry.second)==1,"each face keeps exact source ancestry");
    std::vector<double> lo,hi,xyz;
    gmsh::model::getParametrizationBounds(2,entry.second,lo,hi);
    gmsh::model::getValue(2,entry.second,{(lo[0]+hi[0])/2,(lo[1]+hi[1])/2},xyz);
    require(xyz.size()==3&&std::isfinite(xyz[0]),"each descendant tag is a valid evaluation handle");
  }
}
std::map<std::string, std::string> reader_state() {
  const char *names[] = {
    "read.precision.mode", "read.precision.val", "read.maxprecision.mode",
    "read.maxprecision.val", "read.surfacecurve.mode", "read.stdsameparameter.mode",
    "xstep.cascade.unit", "read.step.shape.repr", "read.step.shape.relationship",
    "read.step.shape.aspect", "read.step.product.mode", "read.step.product.context",
    "read.step.assembly.level", "read.step.nonmanifold", "read.step.ideas",
    "read.step.resource.name", "read.step.sequence", "read.encoderegularity.angle",
    "step.angleunit.mode", "read.step.tessellated", "read.step.constructivegeom.relationship",
    "read.step.codepage", "read.step.all.shapes", "read.step.root.transformation"};
  std::map<std::string, std::string> result;
  for (const char *name : names) {
    const char *value = Interface_Static::CVal(name);
    require(value != nullptr && *value != '\0', std::string("known reader setting missing: ") + name);
    result.emplace(name, value);
  }
  return result;
}
void context_fixture(const std::string &step_file) {
  STEPControl_Controller::Init();
  const auto initial = reader_state();
  stepmesher::occ::Context first(step_file);
  require(reader_state() == initial, "reader state must be restored after default import");
  const auto first_measurements = first.measurements(1e-6);
  const auto original_tags = first.bind_native();
  require(original_tags.size() == first.face_count() && !original_tags.empty(), "complete native face binding");
  require(first.face_occurrences().size() >= first.face_count(), "complete occurrence map");
  std::vector<stepmesher::occ::ProjectNode> points;
  const auto descriptors = first.adapter_face_descriptors();
  for (std::size_t i = 0; i < descriptors.size(); ++i) points.push_back({original_tags[i], descriptors[i].center});
  const auto a = first.project_nodes(points);
  gmsh::clear();
  gmsh::model::add("altered-reader-state");
  require(Interface_Static::SetIVal("read.precision.mode", 1), "set nondefault precision mode");
  require(Interface_Static::SetRVal("read.precision.val", 0.02), "set nondefault precision");
  require(Interface_Static::SetIVal("read.step.ideas", 1), "set ideas import mode");
  require(Interface_Static::SetIVal("read.step.nonmanifold", 1), "set nonmanifold import mode");
  require(Interface_Static::SetCVal("xstep.cascade.unit", "M"), "set nondefault units");
  const auto modified = reader_state();
  stepmesher::occ::Context second(step_file);
  require(reader_state() == modified, "reader state must be restored after altered import");
  const auto second_measurements = second.measurements(1e-6);
  require(first_measurements.diag == second_measurements.diag && first_measurements.area == second_measurements.area,
          "old subprocess defaults must be independent of caller state");
  const auto second_tags = second.bind_native();
  require(original_tags == second_tags, "same native import tags after fresh model");
  const auto b = second.project_nodes(points);
  require(a.n_fallback == b.n_fallback && a.nodes.size() == b.nodes.size(), "projection statuses preserved");
  for (std::size_t i = 0; i < a.nodes.size(); ++i)
    require(a.nodes[i].x == b.nodes[i].x && a.nodes[i].uv == b.nodes[i].uv && a.nodes[i].distance == b.nodes[i].distance,
            "projection coordinates independent of caller reader state");
  for (const auto &entry : initial) Interface_Static::SetCVal(entry.first.c_str(), entry.second.c_str());
  gmsh::clear();
  gmsh::model::add("guarded-options");
  stepmesher::occ::Context guarded(step_file);
  for (const auto &entry : std::vector<std::pair<std::string,double>>{
         {"Geometry.OCCFixDegenerated",1}, {"Geometry.OCCFixSmallEdges",1},
         {"Geometry.OCCFixSmallFaces",1}, {"Geometry.OCCSewFaces",1},
         {"Geometry.OCCMakeSolids",1}, {"Geometry.OCCScaling",2}}) {
    gmsh::option::setNumber(entry.first, entry.second);
    bool rejected = false;
    try { guarded.bind_native(); } catch (const std::exception &e) {
      rejected = std::string(e.what()).find("[native_runtime]") == 0;
    }
    require(rejected, "unsupported import option must fail before passing a pointer: " + entry.first);
    gmsh::option::setNumber(entry.first, entry.first == "Geometry.OCCScaling" ? 1.0 : 0.0);
    require(inventory().empty(), "option rejection must not import geometry");
  }
}
}

void run_occ_context_probe(const std::string &step_file) {
  STEPControl_Controller::Init();
  // These are pinned OCCT 7.9.3 defaults of the old independent helper process.
  require(Interface_Static::IVal("read.precision.mode") == 0, "fresh precision mode");
  require(Interface_Static::RVal("read.precision.val") == 1e-3, "fresh precision value");
  require(Interface_Static::IVal("read.maxprecision.mode") == 0, "fresh maxprecision mode");
  require(Interface_Static::RVal("read.maxprecision.val") == 1, "fresh maxprecision value");
  require(Interface_Static::IVal("read.surfacecurve.mode") == 0, "fresh surfacecurve mode");
  require(Interface_Static::IVal("read.stdsameparameter.mode") == 0, "fresh SameParameter mode");
  require(Interface_Static::IVal("read.step.ideas") == 0 && Interface_Static::IVal("read.step.nonmanifold") == 0,
          "fresh STEP ideas/nonmanifold modes");
  for (int repeat = 0; repeat < 3; ++repeat) {
    gmsh::initialize(0, nullptr, false, false);
    try {
      gmsh::option::setNumber("General.Terminal", 0);
      gmsh::model::add("native-located-shared");
      native_fixture(false);
      gmsh::clear(); gmsh::model::add("native-coincident-descriptors");
      native_fixture(true);
      gmsh::clear(); gmsh::model::add("shared-support-distinct-trims");
      shared_support_trim_fixture();
      gmsh::clear(); gmsh::model::add("periodic-seams");
      periodic_fixture();
      gmsh::clear(); gmsh::model::add("fragment-descendants");
      fragment_fixture();
      gmsh::clear(); gmsh::model::add("retained-context");
      context_fixture(step_file);
      gmsh::finalize();
    } catch (...) {
      gmsh::finalize();
      throw;
    }
  }
}

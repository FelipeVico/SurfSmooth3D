#include "step_mesher/mesher.hpp"

#include "mex.h"

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string mx_to_string(const mxArray *value, const char *name) {
  if (!value || !mxIsChar(value)) {
    throw std::runtime_error(std::string(name) + " must be a char/string");
  }
  char *raw = mxArrayToString(value);
  if (!raw) {
    throw std::runtime_error(std::string("Could not convert ") + name);
  }
  std::string out(raw);
  mxFree(raw);
  return out;
}

double get_scalar_field(const mxArray *s, const char *name, double fallback) {
  if (!s || !mxIsStruct(s)) {
    return fallback;
  }
  const mxArray *field = mxGetField(s, 0, name);
  if (!field || mxIsEmpty(field)) {
    return fallback;
  }
  return mxGetScalar(field);
}

std::string get_string_field(const mxArray *s, const char *name,
                             const std::string &fallback) {
  if (!s || !mxIsStruct(s)) {
    return fallback;
  }
  const mxArray *field = mxGetField(s, 0, name);
  if (!field || mxIsEmpty(field)) {
    return fallback;
  }
  return mx_to_string(field, name);
}

bool get_bool_field(const mxArray *s, const char *name, bool fallback) {
  if (!s || !mxIsStruct(s)) {
    return fallback;
  }
  const mxArray *field = mxGetField(s, 0, name);
  if (!field || mxIsEmpty(field)) {
    return fallback;
  }
  return mxIsLogical(field) ? mxIsLogicalScalarTrue(field)
                            : mxGetScalar(field) != 0.0;
}

stepmesher::MeshOptions parse_options(const mxArray *opts_array) {
  stepmesher::MeshOptions opts;
  opts.order = static_cast<int>(get_scalar_field(opts_array, "order", opts.order));
  opts.mesh_fraction = get_scalar_field(opts_array, "mesh_fraction",
                                        opts.mesh_fraction);
  opts.refinement_level = static_cast<int>(
      get_scalar_field(opts_array, "refinement_level", opts.refinement_level));
  opts.edge_gll_order =
      static_cast<int>(get_scalar_field(opts_array, "edge_gll_order",
                                        opts.edge_gll_order));
  opts.mesh_profile =
      get_string_field(opts_array, "mesh_profile", opts.mesh_profile);
  opts.curvature_points =
      static_cast<int>(get_scalar_field(opts_array, "curvature_points",
                                        opts.curvature_points));
  opts.hmin_mode = get_string_field(opts_array, "hmin_mode", opts.hmin_mode);
  opts.hmin_fraction =
      get_scalar_field(opts_array, "hmin_fraction", opts.hmin_fraction);
  opts.hmin_edge_scale =
      get_scalar_field(opts_array, "hmin_edge_scale", opts.hmin_edge_scale);
  opts.cad_edge_tol_fraction = get_scalar_field(
      opts_array, "cad_edge_tol_fraction", opts.cad_edge_tol_fraction);
  opts.hmax_fraction =
      get_scalar_field(opts_array, "hmax_fraction", opts.hmax_fraction);
  opts.mesh_size_extend_from_boundary = get_bool_field(
      opts_array, "mesh_size_extend_from_boundary",
      opts.mesh_size_extend_from_boundary);
  opts.gmsh_algorithm =
      get_string_field(opts_array, "gmsh_algorithm", opts.gmsh_algorithm);
  opts.optimize = get_bool_field(opts_array, "optimize", opts.optimize);
  opts.min_angle_deg =
      get_scalar_field(opts_array, "min_angle_deg", opts.min_angle_deg);
  opts.max_edge_ratio =
      get_scalar_field(opts_array, "max_edge_ratio", opts.max_edge_ratio);
  opts.enforce_quality =
      get_bool_field(opts_array, "enforce_quality", opts.enforce_quality);
  opts.occt_precision =
      get_scalar_field(opts_array, "occt_precision", opts.occt_precision);
  opts.occt_maxprecision =
      get_scalar_field(opts_array, "occt_maxprecision", opts.occt_maxprecision);
  opts.occ_info_exe =
      get_string_field(opts_array, "occ_info_exe", opts.occ_info_exe);
  opts.sameparameter =
      get_bool_field(opts_array, "sameparameter", opts.sameparameter);
  opts.surfacecurve_mode = get_string_field(opts_array, "surfacecurve_mode",
                                            opts.surfacecurve_mode);
  opts.verbose = get_bool_field(opts_array, "verbose", opts.verbose);
  return opts;
}

stepmesher::SampleNodes parse_sample_nodes(const mxArray *uv_array, int order) {
  if (!uv_array || !mxIsDouble(uv_array) || mxIsComplex(uv_array) ||
      mxGetM(uv_array) != 2) {
    throw std::runtime_error("uv must be a real double array of size 2 x N");
  }

  stepmesher::SampleNodes nodes;
  nodes.order = order;
  nodes.nodes_per_patch = static_cast<int>(mxGetN(uv_array));
  const double *uv = mxGetPr(uv_array);
  nodes.uv.assign(uv, uv + 2 * static_cast<std::size_t>(nodes.nodes_per_patch));
  return nodes;
}

mxArray *make_double_matrix(const std::vector<double> &data, mwSize rows,
                            mwSize cols) {
  mxArray *array = mxCreateDoubleMatrix(rows, cols, mxREAL);
  std::copy(data.begin(), data.end(), mxGetPr(array));
  return array;
}

mxArray *make_double_row_from_ints(const std::vector<int> &data) {
  mxArray *array = mxCreateDoubleMatrix(1, data.size(), mxREAL);
  double *out = mxGetPr(array);
  for (std::size_t i = 0; i < data.size(); ++i) {
    out[i] = static_cast<double>(data[i]);
  }
  return array;
}

mxArray *make_meta_struct(const stepmesher::MeshResult &mesh,
                          const stepmesher::MeshOptions &opts) {
  const char *fields[] = {"cad_identity", "occ_evaluation", "runtime_version",
                          "gmsh_version", "occt_version", "native_identity_faces",
                          "native_identity_occurrences", "backend",
                          "step_file",
                          "mesh_profile",
                          "sample_rule",
                          "mesh_fraction",
                          "refinement_level",
                          "edge_gll_order",
                          "curvature_points",
                          "hmin_mode",
                          "hmin_fraction",
                          "hmin_edge_scale",
                          "cad_edge_tol_fraction",
                          "hmax_fraction",
                          "hmin_effective",
                          "hmax_effective",
                          "mesh_size_extend_from_boundary",
                          "gmsh_algorithm",
                          "optimize",
                          "min_angle_deg",
                          "max_edge_ratio",
                          "enforce_quality",
                          "occt_precision",
                          "occt_maxprecision",
                          "occ_info_exe",
                          "sameparameter",
                          "surfacecurve_mode",
                          "bbox_diag",
                          "mesh_size",
                          "min_raw_cad_edge_length",
                          "min_meaningful_cad_edge_length",
                          "cad_edge_count",
                          "meaningful_cad_edge_count",
                          "cad_area",
                          "max_projection_distance",
                          "mean_projection_distance",
                          "projection_failures",
                          "step1b_edge_newton_success",
                          "step1b_edge_curve_projection_success",
                          "edge_newton_success",
                          "edge_curve_projection_success",
                          "edge_fallback_points",
                          "continuity_break_count",
                          "continuity_helper_curve_count",
                          "fragmented_face_count",
                          "feature_edge_fallback_points",
                          "edge_order",
                          "edge_nodes_per_panel",
                          "ncad_edges",
                          "npanels_total",
                          "min_scaffold_angle_deg",
                          "max_scaffold_edge_ratio",
                          "bad_min_angle_count",
                          "bad_edge_ratio_count",
                          "worst_quality_triangle"};
  mxArray *meta =
      mxCreateStructMatrix(1, 1, sizeof(fields) / sizeof(fields[0]), fields);
  mxSetField(meta, 0, "backend", mxCreateString(mesh.backend.c_str()));
  mxSetField(meta, 0, "cad_identity", mxCreateString(mesh.cad_identity.c_str()));
  mxSetField(meta, 0, "occ_evaluation", mxCreateString(mesh.occ_evaluation.c_str()));
  mxSetField(meta, 0, "runtime_version", mxCreateString(mesh.runtime_version.c_str()));
  mxSetField(meta, 0, "gmsh_version", mxCreateString(mesh.gmsh_version.c_str()));
  mxSetField(meta, 0, "occt_version", mxCreateString(mesh.occt_version.c_str()));
  mxSetField(meta, 0, "native_identity_faces", mxCreateDoubleScalar(mesh.native_identity_faces));
  mxSetField(meta, 0, "native_identity_occurrences", mxCreateDoubleScalar(mesh.native_identity_occurrences));
  mxSetField(meta, 0, "step_file", mxCreateString(mesh.step_file.c_str()));
  mxSetField(meta, 0, "mesh_profile",
             mxCreateString(mesh.mesh_profile.c_str()));
  mxSetField(meta, 0, "sample_rule",
             mxCreateString(mesh.sample_rule.c_str()));
  mxSetField(meta, 0, "mesh_fraction", mxCreateDoubleScalar(opts.mesh_fraction));
  mxSetField(meta, 0, "refinement_level",
             mxCreateDoubleScalar(opts.refinement_level));
  mxSetField(meta, 0, "edge_gll_order",
             mxCreateDoubleScalar(opts.edge_gll_order));
  mxSetField(meta, 0, "curvature_points",
             mxCreateDoubleScalar(opts.curvature_points));
  mxSetField(meta, 0, "hmin_mode", mxCreateString(opts.hmin_mode.c_str()));
  mxSetField(meta, 0, "hmin_fraction",
             mxCreateDoubleScalar(opts.hmin_fraction));
  mxSetField(meta, 0, "hmin_edge_scale",
             mxCreateDoubleScalar(opts.hmin_edge_scale));
  mxSetField(meta, 0, "cad_edge_tol_fraction",
             mxCreateDoubleScalar(opts.cad_edge_tol_fraction));
  mxSetField(meta, 0, "hmax_fraction",
             mxCreateDoubleScalar(opts.hmax_fraction));
  mxSetField(meta, 0, "hmin_effective",
             mxCreateDoubleScalar(mesh.hmin_effective));
  mxSetField(meta, 0, "hmax_effective",
             mxCreateDoubleScalar(mesh.hmax_effective));
  mxSetField(meta, 0, "mesh_size_extend_from_boundary",
             mxCreateLogicalScalar(opts.mesh_size_extend_from_boundary));
  mxSetField(meta, 0, "gmsh_algorithm",
             mxCreateString(opts.gmsh_algorithm.c_str()));
  mxSetField(meta, 0, "optimize", mxCreateLogicalScalar(opts.optimize));
  mxSetField(meta, 0, "min_angle_deg",
             mxCreateDoubleScalar(opts.min_angle_deg));
  mxSetField(meta, 0, "max_edge_ratio",
             mxCreateDoubleScalar(opts.max_edge_ratio));
  mxSetField(meta, 0, "enforce_quality",
             mxCreateLogicalScalar(opts.enforce_quality));
  mxSetField(meta, 0, "occt_precision",
             mxCreateDoubleScalar(opts.occt_precision));
  mxSetField(meta, 0, "occt_maxprecision",
             mxCreateDoubleScalar(opts.occt_maxprecision));
  mxSetField(meta, 0, "occ_info_exe",
             mxCreateString(opts.occ_info_exe.c_str()));
  mxSetField(meta, 0, "sameparameter",
             mxCreateLogicalScalar(opts.sameparameter));
  mxSetField(meta, 0, "surfacecurve_mode",
             mxCreateString(opts.surfacecurve_mode.c_str()));
  mxSetField(meta, 0, "bbox_diag", mxCreateDoubleScalar(mesh.bbox_diag));
  mxSetField(meta, 0, "mesh_size", mxCreateDoubleScalar(mesh.mesh_size));
  mxSetField(meta, 0, "min_raw_cad_edge_length",
             mxCreateDoubleScalar(mesh.min_raw_cad_edge_length));
  mxSetField(meta, 0, "min_meaningful_cad_edge_length",
             mxCreateDoubleScalar(mesh.min_meaningful_cad_edge_length));
  mxSetField(meta, 0, "cad_edge_count",
             mxCreateDoubleScalar(mesh.cad_edge_count));
  mxSetField(meta, 0, "meaningful_cad_edge_count",
             mxCreateDoubleScalar(mesh.meaningful_cad_edge_count));
  mxSetField(meta, 0, "cad_area", mxCreateDoubleScalar(mesh.cad_area));
  mxSetField(meta, 0, "max_projection_distance",
             mxCreateDoubleScalar(mesh.max_projection_distance));
  mxSetField(meta, 0, "mean_projection_distance",
             mxCreateDoubleScalar(mesh.mean_projection_distance));
  mxSetField(meta, 0, "projection_failures",
             mxCreateDoubleScalar(mesh.projection_failures));
  mxSetField(meta, 0, "step1b_edge_newton_success",
             mxCreateDoubleScalar(mesh.step1b_edge_newton_success));
  mxSetField(meta, 0, "step1b_edge_curve_projection_success",
             mxCreateDoubleScalar(mesh.step1b_edge_curve_projection_success));
  mxSetField(meta, 0, "edge_newton_success",
             mxCreateDoubleScalar(mesh.edge_newton_success));
  mxSetField(meta, 0, "edge_curve_projection_success",
             mxCreateDoubleScalar(mesh.edge_curve_projection_success));
  mxSetField(meta, 0, "edge_fallback_points",
             mxCreateDoubleScalar(mesh.edge_fallback_points));
  mxSetField(meta, 0, "continuity_break_count",
             mxCreateDoubleScalar(mesh.continuity_break_count));
  mxSetField(meta, 0, "continuity_helper_curve_count",
             mxCreateDoubleScalar(mesh.continuity_helper_curve_count));
  mxSetField(meta, 0, "fragmented_face_count",
             mxCreateDoubleScalar(mesh.fragmented_face_count));
  mxSetField(meta, 0, "feature_edge_fallback_points",
             mxCreateDoubleScalar(mesh.feature_edge_fallback_points));
  mxSetField(meta, 0, "edge_order", mxCreateDoubleScalar(mesh.edge_order));
  mxSetField(meta, 0, "edge_nodes_per_panel",
             mxCreateDoubleScalar(mesh.edge_nodes_per_panel));
  mxSetField(meta, 0, "ncad_edges", mxCreateDoubleScalar(mesh.ncad_edges));
  mxSetField(meta, 0, "npanels_total",
             mxCreateDoubleScalar(mesh.npanels_total));
  mxSetField(meta, 0, "min_scaffold_angle_deg",
             mxCreateDoubleScalar(mesh.min_scaffold_angle_deg));
  mxSetField(meta, 0, "max_scaffold_edge_ratio",
             mxCreateDoubleScalar(mesh.max_scaffold_edge_ratio));
  mxSetField(meta, 0, "bad_min_angle_count",
             mxCreateDoubleScalar(mesh.bad_min_angle_count));
  mxSetField(meta, 0, "bad_edge_ratio_count",
             mxCreateDoubleScalar(mesh.bad_edge_ratio_count));
  mxSetField(meta, 0, "worst_quality_triangle",
             mxCreateDoubleScalar(mesh.worst_quality_triangle));
  return meta;
}

mxArray *mesh_to_struct(const stepmesher::MeshResult &mesh,
                        const stepmesher::MeshOptions &opts) {
  const char *fields[] = {"xyz",
                          "proj_dist",
                          "nodes_per_patch",
                          "npatches",
                          "norder",
                          "iptype",
                          "tri_face_tags",
                          "parent_coarse_tri",
                          "scaffold_nodes",
                          "scaffold_node_tags",
                          "scaffold_triangles",
                          "scaffold_triangle_node_tags",
                          "scaffold_face_tags",
                          "scaffold_min_angle_deg",
                          "scaffold_edge_ratio",
                          "scaffold_triangle_area",
                          "edge_order",
                          "edge_nodes_per_panel",
                          "ncad_edges",
                          "npanels_total",
                          "edge_gll_nodes",
                          "edge_ids",
                          "edge_curve_tags",
                          "edge_adjacent_faces",
                          "edge_panel_start",
                          "edge_panel_count",
                          "panel_node_tags",
                          "panel_xyz",
                          "panel_tangents",
                          "panel_weights",
                          "meta"};
  mxArray *out =
      mxCreateStructMatrix(1, 1, sizeof(fields) / sizeof(fields[0]), fields);
  mxSetField(out, 0, "xyz",
             make_double_matrix(mesh.xyz, 3,
                                static_cast<mwSize>(mesh.xyz.size() / 3)));
  mxSetField(out, 0, "proj_dist",
             make_double_matrix(mesh.proj_dist, 1,
                                static_cast<mwSize>(mesh.proj_dist.size())));
  mxSetField(out, 0, "nodes_per_patch",
             mxCreateDoubleScalar(mesh.nodes_per_patch));
  mxSetField(out, 0, "npatches", mxCreateDoubleScalar(mesh.npatches));
  mxSetField(out, 0, "norder", mxCreateDoubleScalar(mesh.order));
  mxSetField(out, 0, "iptype", mxCreateDoubleScalar(mesh.iptype));
  mxSetField(out, 0, "tri_face_tags",
             make_double_row_from_ints(mesh.tri_face_tags));
  mxSetField(out, 0, "parent_coarse_tri",
             make_double_row_from_ints(mesh.parent_coarse_tri));
  mxSetField(out, 0, "scaffold_nodes",
             make_double_matrix(mesh.scaffold_nodes, 3,
                                static_cast<mwSize>(mesh.nscaffold_nodes)));
  mxSetField(out, 0, "scaffold_node_tags",
             make_double_row_from_ints(mesh.scaffold_node_tags));
  mxSetField(out, 0, "scaffold_triangles",
             make_double_row_from_ints(mesh.scaffold_triangles));
  mxArray *scaffold_triangles = mxGetField(out, 0, "scaffold_triangles");
  mxSetM(scaffold_triangles, 3);
  mxSetN(scaffold_triangles, static_cast<mwSize>(mesh.nscaffold_triangles));
  mxSetField(out, 0, "scaffold_triangle_node_tags",
             make_double_row_from_ints(mesh.scaffold_triangle_node_tags));
  mxArray *scaffold_triangle_node_tags =
      mxGetField(out, 0, "scaffold_triangle_node_tags");
  mxSetM(scaffold_triangle_node_tags, 3);
  mxSetN(scaffold_triangle_node_tags,
         static_cast<mwSize>(mesh.nscaffold_triangles));
  mxSetField(out, 0, "scaffold_face_tags",
             make_double_row_from_ints(mesh.scaffold_face_tags));
  mxSetField(out, 0, "scaffold_min_angle_deg",
             make_double_matrix(mesh.scaffold_min_angle_deg, 1,
                                static_cast<mwSize>(mesh.scaffold_min_angle_deg.size())));
  mxSetField(out, 0, "scaffold_edge_ratio",
             make_double_matrix(mesh.scaffold_edge_ratio, 1,
                                static_cast<mwSize>(mesh.scaffold_edge_ratio.size())));
  mxSetField(out, 0, "scaffold_triangle_area",
             make_double_matrix(mesh.scaffold_triangle_area, 1,
                                static_cast<mwSize>(mesh.scaffold_triangle_area.size())));
  mxSetField(out, 0, "edge_order", mxCreateDoubleScalar(mesh.edge_order));
  mxSetField(out, 0, "edge_nodes_per_panel",
             mxCreateDoubleScalar(mesh.edge_nodes_per_panel));
  mxSetField(out, 0, "ncad_edges", mxCreateDoubleScalar(mesh.ncad_edges));
  mxSetField(out, 0, "npanels_total",
             mxCreateDoubleScalar(mesh.npanels_total));
  mxSetField(out, 0, "edge_gll_nodes",
             make_double_matrix(mesh.edge_gll_nodes, 1,
                                static_cast<mwSize>(mesh.edge_gll_nodes.size())));
  mxSetField(out, 0, "edge_ids", make_double_row_from_ints(mesh.edge_ids));
  mxSetField(out, 0, "edge_curve_tags",
             make_double_row_from_ints(mesh.edge_curve_tags));
  mxSetField(out, 0, "edge_adjacent_faces",
             make_double_row_from_ints(mesh.edge_adjacent_faces));
  mxArray *edge_adjacent_faces = mxGetField(out, 0, "edge_adjacent_faces");
  mxSetM(edge_adjacent_faces, 2);
  mxSetN(edge_adjacent_faces, static_cast<mwSize>(mesh.ncad_edges));
  mxSetField(out, 0, "edge_panel_start",
             make_double_row_from_ints(mesh.edge_panel_start));
  mxSetField(out, 0, "edge_panel_count",
             make_double_row_from_ints(mesh.edge_panel_count));
  mxSetField(out, 0, "panel_node_tags",
             make_double_row_from_ints(mesh.panel_node_tags));
  mxArray *panel_node_tags = mxGetField(out, 0, "panel_node_tags");
  mxSetM(panel_node_tags, 2);
  mxSetN(panel_node_tags, static_cast<mwSize>(mesh.npanels_total));
  mxSetField(out, 0, "panel_xyz",
             make_double_matrix(mesh.panel_xyz, 3,
                                static_cast<mwSize>(mesh.panel_xyz.size() / 3)));
  mxSetField(out, 0, "panel_tangents",
             make_double_matrix(mesh.panel_tangents, 3,
                                static_cast<mwSize>(mesh.panel_tangents.size() / 3)));
  mxSetField(out, 0, "panel_weights",
             make_double_matrix(mesh.panel_weights, 1,
                                static_cast<mwSize>(mesh.panel_weights.size())));
  mxSetField(out, 0, "meta", make_meta_struct(mesh, opts));
  return out;
}

void mex_mesh_step(int nlhs, mxArray **plhs, int nrhs,
                   const mxArray **prhs) {
  if (nrhs < 4) {
    throw std::runtime_error(
        "mesh_step requires: command, stepFile, opts, uv");
  }
  if (nlhs > 1) {
    throw std::runtime_error("mesh_step returns one output struct");
  }

  const std::string step_file = mx_to_string(prhs[1], "stepFile");
  const stepmesher::MeshOptions opts = parse_options(prhs[2]);
  const stepmesher::SampleNodes samples = parse_sample_nodes(prhs[3], opts.order);
  stepmesher::MeshResult mesh =
      stepmesher::mesh_step_file(step_file, opts, samples);
  mesh.sample_rule = "rv_nodes_from_matlab";
  plhs[0] = mesh_to_struct(mesh, opts);
}

} // namespace

void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs) {
  try {
    if (nrhs < 1) {
      throw std::runtime_error("First argument must be a command string");
    }
    const std::string command = mx_to_string(prhs[0], "command");
    if (command == "mesh_step") {
      mex_mesh_step(nlhs, plhs, nrhs, prhs);
    } else {
      throw std::runtime_error("Unknown step_mesher_mex command: " + command);
    }
  } catch (const std::exception &e) {
    const std::string message = e.what();
    const char *id = "stepmesher:mex:error";
    if (message.find("[native_identity]") != std::string::npos) id = "stepmesher:nativeIdentity";
    else if (message.find("[native_runtime]") != std::string::npos) id = "stepmesher:nativeRuntime";
    else if (message.find("[unsupported_occ_helper]") != std::string::npos) id = "stepmesher:unsupportedOccHelper";
    mexErrMsgIdAndTxt(id, "%s", e.what());
  }
}

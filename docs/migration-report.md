# SurfSmooth3D migration and ownership report

Audit date: 6 October 2026. Scope: the **307 tracked files** in the completed extraction at commit `8b2ff6c0a01410a21296821d3a2f317b3d1c3f0d`. Generated binaries, local dependency installations, build logs and output directories are excluded. This report and its CSV are subsequent documentation, outside that 307-file baseline.

**The package separation is sound, but “every copied file must be unnecessary in fmm3dbie” is too strict.** Geometry generation belongs in SurfSmooth3D. Geometry representation, polynomial evaluation, interpolation and quadrature are needed by both packages. Even a BIE solver that accepts only `.go3` must reconstruct coefficients/weights and evaluate/oversample its input surface.

The extraction **copied** files; it did not delete or modify the original checkouts. Thus there are two different kinds of duplication: temporary old/new copies retained for comparison, and genuinely shared utilities that both packages must continue to have. The extraction also preserved complete tested modules/classes rather than attempting to remove every optional helper. Those exceptions are identified below.

## 1. Complete accounting

The companion `docs/migration-inventory.csv` contains **one row for every baseline file**, with its exact destination, source or generation procedure, purpose, requirement level, fmm3dbie retention decision, duplication explanation and SHA-256. Paths in its `file` column are relative to SurfSmooth3D; original source paths are relative to the containing workspace. A generated gateway source is relative to SurfSmooth3D.

| Origin | Files | What the number means |
|---|---:|---|
| Explicit imports in the sealed extraction manifest | 230 | 162 from fmm3dbie proper, 2 uncompiled reference files from its FMM3D checkout, and 66 from STEP V3. Of these, 116 are byte-identical and 114 have recorded mechanical changes. |
| Additional exact copies | 22 | 5 attribution/license files, 8 CAD fixture files, 7 archived STEP response/map files, and 2 prior validation records. |
| Archived requests with a portable input path | 5 | Only the STEP input path was replaced with a fixture placeholder. |
| Generated sphere fixtures | 3 | Outputs of the existing BIE sphere generator; the generator itself was not imported. |
| Generated MWrap C gateway | 1 | Regenerated from the narrowed smoother interface. |
| New integration files | 46 | 15 packaging/build/docs/setup/example files, 12 test programs/runners and 19 provenance/reference-configuration records. |
| **Total** | **307** | Includes source, tables, fixtures, tests, build files, documentation and evidence; it is not a count of numerical source files. |

The complete functional partition is:

| Block | Files | Why present / ownership |
|---|---:|---|
| Fortran smoothing core | 14 | Generate the smooth surface from CAD quadrature and a launch scaffold; implement sigma, Newton projection and refinement. |
| Shared active Fortran mathematics/support | 7 | Supply polynomial quadrature/transforms, matrix ABI adapters and set/sort operations without linking the BIE library. |
| Uncompiled FMM helper reference copies | 2 | Preserve the exact Legendre and diagnostic-printing source reference used in the original FMM-linked build. |
| MATLAB public interfaces and setup | 10 | Expose the smoother APIs through distinct package/MEX names and provide standalone path discovery. |
| MATLAB runtime compatibility adapters | 3 | Preserve the tested local MEX shutdown workaround without altering Fortran algorithms. |
| Edge-preserving MATLAB workflow | 36 | Load fixed CAD/edge data, match reference coordinates, blend surfaces, diagnose results and export GUI/batch outputs. |
| Adaptive full-smoothing MATLAB workflow | 4 | Drive native adaptive full smoothing, resample the fixed CAD and save/view the result. |
| Adaptive blending MATLAB workflow | 13 | Retain adaptive blending mathematics, independent error checks and native-session management. |
| Shared MATLAB surface class | 19 | Represent, load and manipulate polynomial surface geometry independently of the BIE installation. |
| Shared MATLAB polynomial/interpolation utilities | 26 | Provide RV nodes/weights, basis derivatives/transforms and Legendre check nodes; preserve wider class support. |
| STEP native implementation | 12 | Read CAD, retain native OCCT identity, mesh using Gmsh, sample/project geometry and expose CLI/MEX entry points. |
| STEP MATLAB workflow | 21 | Convert native STEP results to polynomial surfaces, export CAD/edge packages and inspect geometric quality. |
| STEP build, dependency setup and documentation | 6 | Build and install the C++/MATLAB STEP component against a consistent Gmsh/OCCT prefix; document its contract. |
| User drivers and examples | 10 | Provide executable entry points for uniform/adaptive smoothing, blending, STEP meshing and comparisons. |
| Tests and validation runners | 50 | Preserve regressions and verify extraction parity, fixed sources, MATLAB/native interfaces and GO3 interoperability. |
| Geometry and archived-operation fixtures | 38 | Supply portable known CAD, scaffold, sphere and archived STEP responses for examples and regression checks. |
| Standalone build and repository configuration | 7 | Build/install independently using external FMM3D and separate native/MATLAB objects; preserve source bytes and ignore local products. |
| Package documentation and attribution | 6 | Explain installation, architecture, validation and preserved upstream attribution. |
| Provenance and validation evidence | 23 | Record exact origins, working-tree state, hashes, mechanical changes, dependency configuration and test outcomes. |

## 2. Fortran: geometry generation that belongs in SurfSmooth3D

All fourteen files below are in `src/multiscale_mesher/`, copied **byte-for-byte** from the same folder in `tools/fmm3dbie`. No algorithms were rewritten in this extraction. The reviewed gradient fix, one/two-stage behavior, sigma controls and adaptive implementations came from that already-tested source.

| File | Why it belongs here |
|---|---|
| `ModType_Smooth_Surface.f90` | Shared smoother data types for fixed CAD sources, scaffold, launch state and output geometry. |
| `Mod_TreeLRD.f90` | Spatial tree supporting the smoother sigma field; not the BIE near-quadrature tree. |
| `Mod_Fast_Sigma.f90` | Construct/evaluate smoothing length sigma and its gradient, including the reviewed sigma modes. |
| `Mod_Feval.f90` | Evaluate the regularized level-set field and derivatives used by projection. |
| `Mod_Smooth_Surface.f90` | One/two-stage Newton projection, launch-scaffold updates, guards and refinement while preserving fixed CAD sources. |
| `Mod_Adaptive_Smooth_Surface.f90` | Adaptive full smoothing, coefficient-tail/independent checks, leaf refinement and outputs. |
| `Mod_Adaptive_Blend_Surface.f90` | Native session backend for adaptive blending: open/get/refine/project/sigma/close. Blending orchestration remains MATLAB. |
| `surface_smoother.f90` | Public file-based smoother and sigma entry points, filetype dispatch and Fortran C-string adapters. |
| `cisurf_loadmsh.f90` | Load supported scaffold mesh formats and write smoother geometry outputs. |
| `cisurf_skeleton.f90` | Load/build fixed CAD quadrature skeleton and scaffold sources. |
| `cisurf_tritools.f90` | Triangle interpolation/quadrature/refinement utilities used by the smoother, retained as a whole original file. |
| `tfmm_setsub.f` | Smoother FMM adapter, near corrections and regularized small-distance kernel; calls external FMM3D. |
| `Mod_Plot_Tools_sigma.f90` | Legacy sigma/tree plotting diagnostics; kept with the original module group, not called by current workflows. |
| `cisurf_plottools.f90` | Legacy geometry/VTK plotting helpers; retained to preserve the original diagnostic API. |

These files implement creation/refinement of geometry before the GO3 handoff. Static inspection found no calls to their APIs from BIE source outside its `multiscale_mesher` subtree. They are therefore candidates to leave fmm3dbie **together with the corresponding build entries, MATLAB smoother bindings, examples and tests**. This is an ownership conclusion, not an instruction to delete files in isolation.

Two core files are not necessary for the currently exercised production path: `Mod_Plot_Tools_sigma.f90` and `cisurf_plottools.f90` preserve legacy plotting diagnostics. The main source still has `use Mod_Plot_Tools_sigma`, even though no active `plot_sigma` call was found. Removing them would be a separate, small dependency/API cleanup. `cisurf_tritools.f90` was retained whole rather than pruned routine by routine.

## 3. Shared Fortran routines: duplication is currently necessary

| Files in SurfSmooth3D | Why SurfSmooth3D needs them | Why fmm3dbie still needs them |
|---|---|---|
| `src/special_functions/koornexps.f90`, `koorn-uvs-dat.txt`, `koorn-wts-dat.txt` | Koornwinder transforms, derivatives, RV triangle nodes and weights for geometry construction, evaluation and adaptation. | GO3 reading, coefficient reconstruction, interpolation and quadrature. These are mathematical surface tools, not CAD generation. |
| `src/support/lapack_wrap.f90`, `lapack_wrap_64.f90` | Linear algebra with the appropriate native/MATLAB integer ABI. Each build selects one adapter. | The BIE reader and solver infrastructure also use matrix operations. Retain BIE's corresponding adapters. |
| `src/support/setops.f`, `sort.f` | Set operations and sorting used by the smoother FMM/tree interaction. | BIE near-interaction and solver support also use them. |

These **seven source/data files are exact copies**, retained so SurfSmooth3D can build without depending on fmm3dbie. A future special-functions package can own Koornwinder/RV routines; general sorting and BLAS wrappers should receive their own appropriate shared ownership rather than automatically being called “special functions.”

Two more files, `src/special_functions/legeexps.f` and `prini.f`, are **uncompiled reference copies** from the original external FMM source. The makefile deliberately obtains those routines from FMM3D. The copied files are not required for the current build; `prini.f` is printing support, not a special-function algorithm. Their eventual removal or relocation to provenance would not be an algorithm change, but has not been done in this report.

## 4. MATLAB: workflow code and the shared geometry layer

Public MATLAB names are under `surfsmooth3d`. This lets the new workflows and the old BIE MATLAB interfaces coexist without relying on path order.

| Block / files | Purpose | Needed by fmm3dbie after the GO3 boundary? |
|---|---|---|
| `multiscale_mesher.m`, `multiscale_mesher_adaptive.m`, `multiscale_mesher_sigma_eval.m` | One/two-stage smoothing, adaptive full smoothing, and sigma/gradient evaluation. | No, these are geometry-generation APIs. |
| `+edgepreserve/` — 36 files | Fixed CAD/edge input, reference-coordinate matching, blending, sigma controls, validity checks, Newton diagnostics, GUI, saved pairs and batch GO3 export. | No, the finished geometry is consumed by BIE. Plotting/GUI helpers serve the full workflow but are not dependencies of a bare Fortran solve. |
| `+adaptivesmoother/` — 4 files | `solve.m`, `resample_cad.m`, `save_result.m`, `launch_gui.m`. | No. This is adaptive geometry creation, separate from BIE quadrature adaptation. |
| `+adaptiveblend/` — 13 files | Blending basis/beta/CAD evaluation, refinement decisions, independent checks, session orchestration, GUI and export. | No. The blending/adaptation algorithm remains in MATLAB, with a Fortran session backend. |
| `@surfer/` — 19 files | Polynomial surface container, GO3 reader, coefficients/weights, array extraction, plotting and general geometry operations. | **Yes: shared geometry infrastructure.** BIE needs its original class; SurfSmooth3D has a namespaced copy. |
| `+internal/+koorn/` — 6 files | RV nodes/weights and Koornwinder values, derivatives and transforms. | **Yes: shared mathematics.** |
| `+internal/+polytens/` — 19 files | Legendre/Chebyshev polynomial and tensor-product support. | **Yes in BIE's existing broader API.** The main triangular smoother uses a subset; adaptive blending directly needs Legendre check nodes. |
| `+internal/barymat.m` | Barycentric interpolation for the optional Surfacefun adapter. | Shared optional interoperability; not required for ordinary smoothing/GO3 workflows. |

The full class import has explicit exceptions to a strict “only necessary files” rule:

- `@surfer/conv_rsc_to_spmat.m` serves **BIE near-quadrature sparse-matrix assembly**. No smoother workflow/driver caller was found. It belongs in BIE and is a potential later removal from SurfSmooth3D, together with its class declaration.
- `@surfer/surfacemesh_to_surfer.m` and `+internal/barymat.m` preserve optional external Surfacefun/Chebfun interoperability. Ordinary users do not need those external packages.
- `area`, `oversample`, `interpolate_data`, `surf_fun_error`, affine/rotation/scale/translation methods, `merge`, `scatter` and `plot_nodes` preserve the wider class API. They were not all needed or individually exercised by the migrated workflows. In particular BIE still needs its own `oversample`.
- The Chebyshev family and tensor-product parts of the Legendre family preserve quadrilateral support. Do not remove the whole `polytens` folder: `adaptiveblend/basis.m` calls `polytens.lege.pts` even when all output patches are triangular.

Thus the present extraction is a conservative working package, not a proven minimal file set. Optional class pruning should be reviewed separately so the tested numerical workflows stay unchanged.

### Gateways and local runtime adapters

`matlab/surfsmooth3d_adaptive_routs.c` and `surfsmooth3d_adaptive_blend_routs.c` are copies of the dedicated adaptive gateways. `surfsmooth3d_routs.mw` was narrowed to four smoother/filetype signatures, and `surfsmooth3d_routs.c` was generated from it. The broad original `fmm3dbie_routs.c/.mw` **must stay in BIE** because it also exposes solver APIs. Only its smoother-specific bindings can eventually leave.

The three C files in `config/runtime/` preserve the already-tested, platform-conditional MATLAB shutdown workaround. Their purpose is runtime compatibility; they do not implement new smoothing mathematics. Whether the old BIE installation still needs its own workaround is a separate platform question.

`matlab/setup_surfsmooth3d.m`, `+surfsmooth3d/root.m` and the top-level `Contents.m` are new package setup/help files. MWrap itself is needed only to regenerate interfaces; generated C is supplied.

## 5. STEP V3: a separate source project, now part of geometry generation

STEP code came from the **actual working tree** of `step_mesher_cpp_matlab_astra_v3`, including its recorded uncommitted source. It was not extracted from the BIE solver. Its C++ numerical implementation was preserved; the changed C++ text concerns the owned-helper installation path.

| Block | Files / purpose | BIE requirement |
|---|---|---|
| Native implementation — 12 files | `cpp/src/{mesher.cpp,occ_context.cpp,occ_info.cpp,occ_operation.hpp,sample_nodes.cpp,cli.cpp}`, five public headers under `cpp/include/step_mesher/`, and `mex/step_mesher_mex.cpp`. Import CAD, retain topology/identity, mesh/sample/project, and expose native/MATLAB entry points. | None for GO3-only consumption. |
| MATLAB STEP API — 21 files under `matlab/+surfsmooth3d/+stepmesher/` | Run STEP meshing; convert to surface arrays/class; export GO3, scaffold, point quadrature and GLL edges; load/view packages; inspect projection/quality/area/flux; preserve profiles/defaults. | None for the BIE solver. |
| Build/bootstrap/docs — 6 files | `CMakeLists.txt`, `cmake/FindGmsh.cmake`, `dependencies.lock.json`, `tools/build_native_dependencies.py`, `README.md`, `docs/V3_NATIVE_IDENTITY_MIGRATION.md`. | CAD build dependencies belong to SurfSmooth3D. |
| STEP drivers — 4 files | `run_step_mesher.m`, `run_convergence_study.m`, `driver_view_geom_package.m`, `driver_view_step_file.m`. | Geometry-generation examples, not solver dependencies. |
| Tests and fixtures | Nine native test/probe files, three MATLAB test drivers, two CMake check runners, twelve STEP models and thirteen archived-operation fixture/provenance files. | Validation support only. |

The meshing runtime is C++/MATLAB and uses **external Gmsh and OCCT**. The Python file is an optional dependency-build helper, not a runtime mesher. The owned `occ_info` helper is retained: its presence/canonical path participates in selection of the native CAD path, so lack of a subprocess call does not make it safely removable.

CAD edge quadrature, original source surfaces/scaffolds, metadata and parent/reference mappings are needed before and during smoothing/blending. The final BIE GO3 consumer does not need those sidecars for the tested handoff. GO3 alone does not contain connectivity or material information required by every possible BIE formulation.

## 6. Drivers, tests, data and package support

The ten user examples consist of the four STEP drivers above, a new native `examples/fortran/smooth_surface.f90`, and these five migrated MATLAB drivers:

- `driver_two_stage_edge_smooth_batch.m`: two-stage uniform full/partial GUI and batch exports.
- `driver_two_stage_adaptive_smooth.m`: adaptive full smoothing.
- `driver_two_stage_adaptive_blend.m`: adaptive partial blending.
- `driver_edge_preserve_compare_surfaces.m`: comparison/inspection of surfaces.
- `driver_debug_primary_vs_cad_smooth.m`: scaffold-versus-fixed-CAD smoothing diagnostics.

There are **50 test/program/runner files** across native smoother tests, MATLAB tests, interop tests, STEP native probes, STEP MATLAB tests and CMake runners. These belong with the behavior they validate. They are not linked production requirements. The seven `tests/interop/` files combine standalone installed-pipeline checks, result-file comparisons and external BIE consumer/reference/coexistence checks. Where used, the old package is an **external test dependency**, not a runtime dependency of SurfSmooth3D. `test_step_to_smoother.m` explicitly excludes the old BIE implementations from its MATLAB path.

The **38 fixture files** consist of two cylinder inputs, eight CAD package files at two resolutions, three generated sphere GO3 files, twelve STEP models and thirteen archived-operation data/provenance files. The two CAD fixture resolutions intentionally repeat the same scaffold and edge data. This preserves self-contained fixture packages; it is data duplication, not duplicate algorithms. The BIE `+geometries/sphere.m` generator was used to create the sphere fixtures but was not copied into the package.

Seven standalone build/repository configuration files and six package documentation/attribution files make the new package usable independently. The 23 files under `provenance/` preserve hashes, source maps, snapshots, mechanical diffs, attribution and validation evidence. Their presence is justified by reproducibility, not by numerical execution. The CSV distinguishes newly written integration files from copied source, generated files and copied historical evidence.

FMM3D remains an **external numerical dependency**. Neither the BIE solvers nor their layer-potential kernels, singular quadrature tables or GMRES implementations were migrated into SurfSmooth3D. The STEP dependency source trees and compiled binaries are also outside the tracked package baseline.

## 7. What must remain available to fmm3dbie

The ownership boundary is: **SurfSmooth3D constructs geometry; fmm3dbie evaluates and integrates supplied geometry.** fmm3dbie still needs:

1. Its GO3 readers and surface representation.
2. Basis evaluation, differentiation, coefficient conversion, quadrature weights, normals/metrics and interpolation.
3. Oversampling, near-interaction searches, special quadrature and any adaptation needed by its solvers.
4. Shared matrix/set/sort utilities and FMM3D.
5. Its solver gateways, integral kernels, linear solvers and relevant tests.

Concrete existing call sites establish this boundary:

- The Fortran GO3 reader calls `vioreanu_simplex_quad` at line 98, `get_qwts` at line 100 and `dmatvec` at line 106. [Source](/Users/fevibon/Library/CloudStorage/Dropbox/felipe_code_cloud_dropbox/chunkiefest2026/tools/fmm3dbie/src/surface_routs/in_go3.f90:98).
- The MATLAB surface constructor obtains RV nodes/weights and polynomial transforms. [Source](/Users/fevibon/Library/CloudStorage/Dropbox/felipe_code_cloud_dropbox/chunkiefest2026/tools/fmm3dbie/matlab/@surfer/surfer.m:184).
- The Helmholtz solver oversamples the supplied surface. [Source](/Users/fevibon/Library/CloudStorage/Dropbox/felipe_code_cloud_dropbox/chunkiefest2026/tools/fmm3dbie/matlab/+helm3d/+dirichlet/solver.m:116).
- The BIE quadrature-correction builder calls `conv_rsc_to_spmat`. [Source](/Users/fevibon/Library/CloudStorage/Dropbox/felipe_code_cloud_dropbox/chunkiefest2026/tools/fmm3dbie/matlab/+helm3d/+dirichlet/get_quadrature_correction.m:189).
- Adaptive blending in SurfSmooth3D calls Legendre check-node generation. [Source](/Users/fevibon/Library/CloudStorage/Dropbox/felipe_code_cloud_dropbox/chunkiefest2026/SurfSmooth3D/matlab/+surfsmooth3d/+adaptiveblend/basis.m:10).

The future shared special-functions package can remove much of the mathematical source duplication. A shared geometry class is a separate API/ownership decision. For now, separate MATLAB namespaces preserve independence; unchanged native Fortran symbols mean co-linking both native libraries into one executable is not a supported contract of this first extraction.

## 8. Audit result and later cleanup candidates

The main CAD, smoothing, two-stage, edge-preserving, sigma and adaptive workflows are in the correct package. **It would be inaccurate to say all 307 files are necessary for runtime or all can disappear from fmm3dbie.** The CSV marks the distinctions explicitly.

Candidates for a later, separately tested cleanup are the BIE-specific sparse-conversion method in the copied class, optional Surfacefun and wider geometry-class methods, conditional quadrilateral basis support, the two unused plotting modules, and the two uncompiled FMM helper source copies. No removal is proposed for the active shared mathematics without a replacement dependency.

This audit changes documentation only. Counts and origins were reconciled against the committed tree and manifests; extra copied notices/evidence/fixtures were checked against their source bytes. The 230-file extraction hash check was rerun. Numerical algorithms and both original repositories were left unchanged; numerical tests were not repeated solely for this documentation update.

## Appendix: complete file membership by block

Every baseline file appears once below. Use `migration-inventory.csv` for its individual origin, requirement, BIE retention and duplication columns. Names in each block are relative to SurfSmooth3D.

### Fortran smoothing core (14)

- `src/multiscale_mesher/ModType_Smooth_Surface.f90`
- `src/multiscale_mesher/Mod_Adaptive_Blend_Surface.f90`
- `src/multiscale_mesher/Mod_Adaptive_Smooth_Surface.f90`
- `src/multiscale_mesher/Mod_Fast_Sigma.f90`
- `src/multiscale_mesher/Mod_Feval.f90`
- `src/multiscale_mesher/Mod_Plot_Tools_sigma.f90`
- `src/multiscale_mesher/Mod_Smooth_Surface.f90`
- `src/multiscale_mesher/Mod_TreeLRD.f90`
- `src/multiscale_mesher/cisurf_loadmsh.f90`
- `src/multiscale_mesher/cisurf_plottools.f90`
- `src/multiscale_mesher/cisurf_skeleton.f90`
- `src/multiscale_mesher/cisurf_tritools.f90`
- `src/multiscale_mesher/surface_smoother.f90`
- `src/multiscale_mesher/tfmm_setsub.f`

### Shared active Fortran mathematics/support (7)

- `src/special_functions/koorn-uvs-dat.txt`
- `src/special_functions/koorn-wts-dat.txt`
- `src/special_functions/koornexps.f90`
- `src/support/lapack_wrap.f90`
- `src/support/lapack_wrap_64.f90`
- `src/support/setops.f`
- `src/support/sort.f`

### Uncompiled FMM helper reference copies (2)

- `src/special_functions/legeexps.f`
- `src/special_functions/prini.f`

### MATLAB public interfaces and setup (10)

- `matlab/+surfsmooth3d/Contents.m`
- `matlab/+surfsmooth3d/multiscale_mesher.m`
- `matlab/+surfsmooth3d/multiscale_mesher_adaptive.m`
- `matlab/+surfsmooth3d/multiscale_mesher_sigma_eval.m`
- `matlab/+surfsmooth3d/root.m`
- `matlab/setup_surfsmooth3d.m`
- `matlab/surfsmooth3d_adaptive_blend_routs.c`
- `matlab/surfsmooth3d_adaptive_routs.c`
- `matlab/surfsmooth3d_routs.c`
- `matlab/surfsmooth3d_routs.mw`

### MATLAB runtime compatibility adapters (3)

- `config/runtime/surfsmooth3d_adaptive_blend_routs_gateway.c`
- `config/runtime/surfsmooth3d_adaptive_routs_gateway.c`
- `config/runtime/surfsmooth3d_routs_gateway.c`

### Edge-preserving MATLAB workflow (36)

- `matlab/+surfsmooth3d/+edgepreserve/Contents.m`
- `matlab/+surfsmooth3d/+edgepreserve/bad_nodes_table.m`
- `matlab/+surfsmooth3d/+edgepreserve/blend.m`
- `matlab/+surfsmooth3d/+edgepreserve/capture_newton_failure.m`
- `matlab/+surfsmooth3d/+edgepreserve/color_limits.m`
- `matlab/+surfsmooth3d/+edgepreserve/color_values.m`
- `matlab/+surfsmooth3d/+edgepreserve/compute_pair.m`
- `matlab/+surfsmooth3d/+edgepreserve/default_input_root.m`
- `matlab/+surfsmooth3d/+edgepreserve/default_settings.m`
- `matlab/+surfsmooth3d/+edgepreserve/export_surface.m`
- `matlab/+surfsmooth3d/+edgepreserve/file_sha256.m`
- `matlab/+surfsmooth3d/+edgepreserve/launch_gui.m`
- `matlab/+surfsmooth3d/+edgepreserve/load_edges.m`
- `matlab/+surfsmooth3d/+edgepreserve/load_package.m`
- `matlab/+surfsmooth3d/+edgepreserve/mesher_options.m`
- `matlab/+surfsmooth3d/+edgepreserve/new_run_directory.m`
- `matlab/+surfsmooth3d/+edgepreserve/parse_integer_list.m`
- `matlab/+surfsmooth3d/+edgepreserve/plot_gidmsh.m`
- `matlab/+surfsmooth3d/+edgepreserve/plot_newton_failure.m`
- `matlab/+surfsmooth3d/+edgepreserve/plot_surface.m`
- `matlab/+surfsmooth3d/+edgepreserve/print_validity_report.m`
- `matlab/+surfsmooth3d/+edgepreserve/read_gidmsh.m`
- `matlab/+surfsmooth3d/+edgepreserve/read_newton_failure.m`
- `matlab/+surfsmooth3d/+edgepreserve/reference_map.m`
- `matlab/+surfsmooth3d/+edgepreserve/resample_cad.m`
- `matlab/+surfsmooth3d/+edgepreserve/run_batch.m`
- `matlab/+surfsmooth3d/+edgepreserve/save_pair.m`
- `matlab/+surfsmooth3d/+edgepreserve/sigma_mode_label.m`
- `matlab/+surfsmooth3d/+edgepreserve/srcvals_from_positions.m`
- `matlab/+surfsmooth3d/+edgepreserve/stage_source.m`
- `matlab/+surfsmooth3d/+edgepreserve/surface_filename.m`
- `matlab/+surfsmooth3d/+edgepreserve/triangle_maps.m`
- `matlab/+surfsmooth3d/+edgepreserve/update_blend.m`
- `matlab/+surfsmooth3d/+edgepreserve/validity_report.m`
- `matlab/+surfsmooth3d/+edgepreserve/write_go3.m`
- `matlab/+surfsmooth3d/+edgepreserve/write_point_source.m`

### Adaptive full-smoothing MATLAB workflow (4)

- `matlab/+surfsmooth3d/+adaptivesmoother/launch_gui.m`
- `matlab/+surfsmooth3d/+adaptivesmoother/resample_cad.m`
- `matlab/+surfsmooth3d/+adaptivesmoother/save_result.m`
- `matlab/+surfsmooth3d/+adaptivesmoother/solve.m`

### Adaptive blending MATLAB workflow (13)

- `matlab/+surfsmooth3d/+adaptiveblend/Contents.m`
- `matlab/+surfsmooth3d/+adaptiveblend/assess_checks.m`
- `matlab/+surfsmooth3d/+adaptiveblend/basis.m`
- `matlab/+surfsmooth3d/+adaptiveblend/beta_values.m`
- `matlab/+surfsmooth3d/+adaptiveblend/cad_values.m`
- `matlab/+surfsmooth3d/+adaptiveblend/check_error.m`
- `matlab/+surfsmooth3d/+adaptiveblend/default_settings.m`
- `matlab/+surfsmooth3d/+adaptiveblend/independent_indicators.m`
- `matlab/+surfsmooth3d/+adaptiveblend/launch_gui.m`
- `matlab/+surfsmooth3d/+adaptiveblend/nodal_patch.m`
- `matlab/+surfsmooth3d/+adaptiveblend/save_result.m`
- `matlab/+surfsmooth3d/+adaptiveblend/solve.m`
- `matlab/+surfsmooth3d/+adaptiveblend/validate_settings.m`

### Shared MATLAB surface class (19)

- `matlab/+surfsmooth3d/@surfer/Contents.m`
- `matlab/+surfsmooth3d/@surfer/affine_transf.m`
- `matlab/+surfsmooth3d/@surfer/area.m`
- `matlab/+surfsmooth3d/@surfer/conv_rsc_to_spmat.m`
- `matlab/+surfsmooth3d/@surfer/extract_arrays.m`
- `matlab/+surfsmooth3d/@surfer/interpolate_data.m`
- `matlab/+surfsmooth3d/@surfer/load_from_file.m`
- `matlab/+surfsmooth3d/@surfer/merge.m`
- `matlab/+surfsmooth3d/@surfer/oversample.m`
- `matlab/+surfsmooth3d/@surfer/plot.m`
- `matlab/+surfsmooth3d/@surfer/plot_nodes.m`
- `matlab/+surfsmooth3d/@surfer/rotate.m`
- `matlab/+surfsmooth3d/@surfer/scale.m`
- `matlab/+surfsmooth3d/@surfer/scatter.m`
- `matlab/+surfsmooth3d/@surfer/surf_fun_error.m`
- `matlab/+surfsmooth3d/@surfer/surfacemesh_to_surfer.m`
- `matlab/+surfsmooth3d/@surfer/surfer.m`
- `matlab/+surfsmooth3d/@surfer/translate.m`
- `matlab/+surfsmooth3d/@surfer/vals2coefs.m`

### Shared MATLAB polynomial/interpolation utilities (26)

- `matlab/+surfsmooth3d/+internal/+koorn/coefs2vals.m`
- `matlab/+surfsmooth3d/+internal/+koorn/ders.m`
- `matlab/+surfsmooth3d/+internal/+koorn/pols.m`
- `matlab/+surfsmooth3d/+internal/+koorn/rv_nodes.m`
- `matlab/+surfsmooth3d/+internal/+koorn/rv_weights.m`
- `matlab/+surfsmooth3d/+internal/+koorn/vals2coefs.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+cheb/coefs2vals.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+cheb/ders.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+cheb/nodes.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+cheb/pols.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+cheb/pts.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+cheb/vals2coefs.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+cheb/weights.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/coefs2vals.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/ders.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/nodes.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/pol.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/pols.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/polsum.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/pts.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/rts.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/rts_stab.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/tayl.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/vals2coefs.m`
- `matlab/+surfsmooth3d/+internal/+polytens/+lege/weights.m`
- `matlab/+surfsmooth3d/+internal/barymat.m`

### STEP native implementation (12)

- `step_mesher/cpp/include/step_mesher/mesh_result.hpp`
- `step_mesher/cpp/include/step_mesher/mesher.hpp`
- `step_mesher/cpp/include/step_mesher/occ_context.hpp`
- `step_mesher/cpp/include/step_mesher/options.hpp`
- `step_mesher/cpp/include/step_mesher/sample_nodes.hpp`
- `step_mesher/cpp/src/cli.cpp`
- `step_mesher/cpp/src/mesher.cpp`
- `step_mesher/cpp/src/occ_context.cpp`
- `step_mesher/cpp/src/occ_info.cpp`
- `step_mesher/cpp/src/occ_operation.hpp`
- `step_mesher/cpp/src/sample_nodes.cpp`
- `step_mesher/mex/step_mesher_mex.cpp`

### STEP MATLAB workflow (21)

- `matlab/+surfsmooth3d/+stepmesher/apply_mesh_profile_choice.m`
- `matlab/+surfsmooth3d/+stepmesher/export_edges_gll.m`
- `matlab/+surfsmooth3d/+stepmesher/export_geom_package.m`
- `matlab/+surfsmooth3d/+stepmesher/export_gidmsh.m`
- `matlab/+surfsmooth3d/+stepmesher/export_go3.m`
- `matlab/+surfsmooth3d/+stepmesher/export_xyzw.m`
- `matlab/+surfsmooth3d/+stepmesher/load_edges_gll.m`
- `matlab/+surfsmooth3d/+stepmesher/load_geom_package.m`
- `matlab/+surfsmooth3d/+stepmesher/mesh_profile_choices.m`
- `matlab/+surfsmooth3d/+stepmesher/mesh_step.m`
- `matlab/+surfsmooth3d/+stepmesher/plot_coarse_triangle_error.m`
- `matlab/+surfsmooth3d/+stepmesher/plot_geom_package.m`
- `matlab/+surfsmooth3d/+stepmesher/plot_projected_nodes.m`
- `matlab/+surfsmooth3d/+stepmesher/plot_projection_error.m`
- `matlab/+surfsmooth3d/+stepmesher/plot_scaffold.m`
- `matlab/+surfsmooth3d/+stepmesher/plot_scaffold_quality.m`
- `matlab/+surfsmooth3d/+stepmesher/plot_surface.m`
- `matlab/+surfsmooth3d/+stepmesher/private/default_options.m`
- `matlab/+surfsmooth3d/+stepmesher/to_srcvals.m`
- `matlab/+surfsmooth3d/+stepmesher/to_surfer.m`
- `matlab/+surfsmooth3d/+stepmesher/validate_area_flux.m`

### STEP build, dependency setup and documentation (6)

- `step_mesher/CMakeLists.txt`
- `step_mesher/README.md`
- `step_mesher/cmake/FindGmsh.cmake`
- `step_mesher/dependencies.lock.json`
- `step_mesher/docs/V3_NATIVE_IDENTITY_MIGRATION.md`
- `step_mesher/tools/build_native_dependencies.py`

### User drivers and examples (10)

- `examples/fortran/smooth_surface.f90`
- `examples/matlab/driver_debug_primary_vs_cad_smooth.m`
- `examples/matlab/driver_edge_preserve_compare_surfaces.m`
- `examples/matlab/driver_two_stage_adaptive_blend.m`
- `examples/matlab/driver_two_stage_adaptive_smooth.m`
- `examples/matlab/driver_two_stage_edge_smooth_batch.m`
- `examples/step/driver_view_geom_package.m`
- `examples/step/driver_view_step_file.m`
- `examples/step/run_convergence_study.m`
- `examples/step/run_step_mesher.m`

### Tests and validation runners (50)

- `step_mesher/cmake/CheckFrozenOperation.cmake`
- `step_mesher/cmake/CheckRetention.cmake`
- `step_mesher/cpp/tests/native_contract.cpp`
- `step_mesher/cpp/tests/occ_context_probe.cpp`
- `step_mesher/cpp/tests/occ_context_probe.hpp`
- `step_mesher/cpp/tests/occ_context_probe_mex.cpp`
- `step_mesher/cpp/tests/smoke.cpp`
- `step_mesher/cpp/tests/test_occ_context.cpp`
- `step_mesher/cpp/tests/test_occ_frozen.cpp`
- `step_mesher/cpp/tests/test_occ_import_equivalence.cpp`
- `step_mesher/cpp/tests/test_occ_retention.cpp`
- `tests/fortran/test_adaptive_blend.f90`
- `tests/fortran/test_adaptive_smoother.f90`
- `tests/fortran/test_fixed_sources.f90`
- `tests/fortran/test_newton_radius_guard.f90`
- `tests/fortran/test_sigma_modes.f90`
- `tests/fortran/test_surfsmooth.c`
- `tests/fortran/test_surfsmooth.f90`
- `tests/interop/bie_read_go3.f90`
- `tests/interop/capture_equivalence.m`
- `tests/interop/compare_equivalence.m`
- `tests/interop/test_bie_consumer.m`
- `tests/interop/test_bie_readers.m`
- `tests/interop/test_gateway_coexistence.m`
- `tests/interop/test_step_to_smoother.m`
- `tests/matlab/adaptive_blend_lock_baseline.m`
- `tests/matlab/check_adaptive_code.m`
- `tests/matlab/fixture_sphere.m`
- `tests/matlab/run_surfsmooth3d_tests.m`
- `tests/matlab/smoke_adaptive_blend_gui.m`
- `tests/matlab/smoke_adaptive_gui.m`
- `tests/matlab/smoke_edgepreserve_gui_batch.m`
- `tests/matlab/smoke_edgepreserve_real_package.m`
- `tests/matlab/smoke_newton_projection_failure_gui.m`
- `tests/matlab/smoke_newton_radius_guard_gui.m`
- `tests/matlab/smoke_sigma_modes_gui.m`
- `tests/matlab/test_adaptive_accuracy.m`
- `tests/matlab/test_adaptive_blend.m`
- `tests/matlab/test_adaptive_blend_accuracy.m`
- `tests/matlab/test_adaptive_blend_math.m`
- `tests/matlab/test_adaptive_reference_maps.m`
- `tests/matlab/test_adaptive_smoother.m`
- `tests/matlab/test_edgepreserve_fixed_cad.m`
- `tests/matlab/test_newton_guard_batch.m`
- `tests/matlab/test_newton_projection_failure.m`
- `tests/matlab/test_newton_radius_guard.m`
- `tests/matlab/test_sigma_modes.m`
- `tests/step/regression_sweep_examples.m`
- `tests/step/smoke_mesh_step.m`
- `tests/step/test_extraction_parity.m`

### Geometry and archived-operation fixtures (38)

- `step_mesher/examples/step_files/Multi_res.step`
- `step_mesher/examples/step_files/Multiscale_step.step`
- `step_mesher/examples/step_files/UV_curved_bezier.step`
- `step_mesher/examples/step_files/cylinder_bezier_v2.step`
- `step_mesher/examples/step_files/intersection_two_balls.step`
- `step_mesher/examples/step_files/piecewise_bezier_cap.step`
- `step_mesher/examples/step_files/spheres_intersect.step`
- `step_mesher/examples/step_files/test_curved_uv.step`
- `step_mesher/examples/step_files/test_geometry_1.step`
- `step_mesher/examples/step_files/test_geometry_2.step`
- `step_mesher/examples/step_files/test_geometry_manas.step`
- `step_mesher/examples/step_files/torus_ellipse.step`
- `tests/fixtures/cad/convergence_intersection_two_balls/p4_mf0p10000_rl0/edges_gll.txt`
- `tests/fixtures/cad/convergence_intersection_two_balls/p4_mf0p10000_rl0/metadata.txt`
- `tests/fixtures/cad/convergence_intersection_two_balls/p4_mf0p10000_rl0/scaffold.gidmsh`
- `tests/fixtures/cad/convergence_intersection_two_balls/p4_mf0p10000_rl0/surface.go3`
- `tests/fixtures/cad/convergence_intersection_two_balls/p8_mf0p10000_rl1/edges_gll.txt`
- `tests/fixtures/cad/convergence_intersection_two_balls/p8_mf0p10000_rl1/metadata.txt`
- `tests/fixtures/cad/convergence_intersection_two_balls/p8_mf0p10000_rl1/scaffold.gidmsh`
- `tests/fixtures/cad/convergence_intersection_two_balls/p8_mf0p10000_rl1/surface.go3`
- `tests/fixtures/cylinder/cylinder.gidmsh`
- `tests/fixtures/cylinder/cylinder_skeleton.txt`
- `tests/fixtures/spheres/sphere_offset_na1_p4.go3`
- `tests/fixtures/spheres/sphere_unit_na2_p10.go3`
- `tests/fixtures/spheres/sphere_unit_na2_p8.go3`
- `tests/step/fixtures/occ_operations/cad-measurements.expected.txt`
- `tests/step/fixtures/occ_operations/edge-samples.expected.txt`
- `tests/step/fixtures/occ_operations/edge-samples.request.in`
- `tests/step/fixtures/occ_operations/feature-nodes.expected.txt`
- `tests/step/fixtures/occ_operations/feature-nodes.request.in`
- `tests/step/fixtures/occ_operations/feature-samples.expected.txt`
- `tests/step/fixtures/occ_operations/feature-samples.request.in`
- `tests/step/fixtures/occ_operations/frozen-map.txt`
- `tests/step/fixtures/occ_operations/helper-curves.expected.txt`
- `tests/step/fixtures/occ_operations/helper-curves.request.in`
- `tests/step/fixtures/occ_operations/project-nodes.expected.txt`
- `tests/step/fixtures/occ_operations/project-nodes.request.in`
- `tests/step/fixtures/occ_operations/provenance.json`

### Standalone build and repository configuration (7)

- `.gitattributes`
- `.gitignore`
- `config/README.md`
- `config/apple-silicon.mk`
- `config/linux.mk`
- `config/runtime/README.md`
- `makefile`

### Package documentation and attribution (6)

- `COPYRIGHT`
- `LICENSE`
- `README.md`
- `THIRD_PARTY.md`
- `docs/architecture.md`
- `docs/validation.md`

### Provenance and validation evidence (23)

- `provenance/build-environment.json`
- `provenance/coexistence-validation.json`
- `provenance/extraction-manifest.json`
- `provenance/fmm-reference-build.json`
- `provenance/fmm-reference-make.inc`
- `provenance/fortran-imports.json`
- `provenance/integration-validation.json`
- `provenance/licenses/FMM3D-AUTHORS`
- `provenance/licenses/FMM3D-LICENSE`
- `provenance/licenses/fmm3dbie-AUTHORS`
- `provenance/matlab-fixtures.json`
- `provenance/matlab-imports.json`
- `provenance/matlab-validation.json`
- `provenance/mechanical-changes.patch`
- `provenance/mwrap-reproducibility.json`
- `provenance/native-validation.json`
- `provenance/prior-migration-legacy-comparison.json`
- `provenance/prior-migration-summary.json`
- `provenance/source-snapshot.json`
- `provenance/step-imports.json`
- `provenance/step-validation.json`
- `provenance/step_mesher-working-tree.patch`
- `provenance/verify_sources.py`


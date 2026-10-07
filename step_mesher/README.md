# SurfSmooth3D STEP component

This optional component imports STEP and constructs the CAD package consumed by
SurfSmooth3D. It copies the tested retained-OCC V3 implementation, including its
working-tree changes. Its meshing, projection, fallback, GLL, blending and manual
refinement algorithms are unchanged. The C++ component version remains `3.0.0`;
its shared library is named `surfsmooth3d_step_core` to avoid an old V3 library
being selected accidentally.

The smoother can consume existing CAD packages without this component. STEP
support needs Gmsh 4.12.2 and OCCT 7.9.3 built together in a common prefix. A
separately installed Gmsh or a second OCCT copy is not a compatible substitute for
native topology sharing. The mesher itself has no Python runtime dependency.
The optional dependency bootstrap uses Python 3.9+, CMake 3.20+, a C++17 compiler,
make, and Freetype development files for Gmsh CAF detection.

From the SurfSmooth3D repository root, build the pinned dependencies offline:

```sh
python3 step_mesher/tools/build_native_dependencies.py \
  --build-root build/deps/cad --archives-dir /path/to/verified/archives --jobs 6
```

`step_mesher/dependencies.lock.json` records the archive names, hashes, versions
and licenses. Supply `--download` instead of `--archives-dir` to fetch missing
archives. The bootstrap never silently downloads dependencies. Native wrappers
are disabled in Gmsh's dependency build; MATLAB uses this component's own MEX.
The supplied recipe is tested on macOS Apple Silicon. Linux verification is not
claimed.

Build and install the component, also from the repository root:

```sh
cmake -S step_mesher -B build/step \
  -DSTEP_MESHER_DEPENDENCY_PREFIX="$PWD/build/deps/cad/install" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD" \
  -DMatlab_ROOT_DIR=/Applications/MATLAB_R2026a.app
cmake --build build/step --parallel 4
ctest --test-dir build/step --output-on-failure
cmake --install build/step
```

Use `STEP_MESHER_BUILD_MEX=OFF` for a native-only build. Use
`STEP_MESHER_SANITIZERS=ON` for ASan/UBSan on the application, not on the external
Gmsh/OCCT libraries. Configuration and compilation stage outputs under
`build/step`; installation publishes the library under `lib`, the diagnostic
CLI under `bin`, and MEX/helper under
`matlab/+surfsmooth3d/+stepmesher/private`. An arbitrary installation prefix is
supported, while the compatible dependency prefix must stay at its configured
location. Publish the complete package, not just its MEX. Start a fresh MATLAB
process after replacing a compiled package.

## MATLAB and CAD packages

```matlab
addpath('/path/to/SurfSmooth3D/matlab');
setup_surfsmooth3d;
opts = struct('order',8, 'edge_gll_order',8, ...
    'mesh_fraction',0.10, 'refinement_level',0, 'mesh_profile','rigid');
mesh = surfsmooth3d.stepmesher.mesh_step( ...
    'step_mesher/examples/step_files/test_curved_uv.step',opts);
[srcvals,norders,iptype,meta] = surfsmooth3d.stepmesher.to_srcvals(mesh);
stats = surfsmooth3d.stepmesher.validate_area_flux(mesh);
pkg = surfsmooth3d.stepmesher.export_geom_package(mesh, ...
    'outputs/step/test_curved_uv');
```

No fmm3dbie MATLAB path is needed. The copied RV/polynomial routines are under
`surfsmooth3d.internal.koorn`; the surface class is `surfsmooth3d.surfer`.
`surface.go3`, `scaffold.gidmsh`, `edges_gll.txt`, and `metadata.txt` retain the
original formats and the reference convention
`X(u,v)=(1-u-v)*P1+u*P2+v*P3`. GO3 uses the polynomial degree directly in its
first line, the number of triangular patches next, then the twelve
component-major blocks of positions, local tangents, and normals at RV nodes.

The ordinary driver is `examples/step/run_step_mesher.m`; it currently uses
p4/GLL4 and `curvature_balanced_local`. Uncomment a `stepFile` selection in the
driver before running it. API defaults remain p8/GLL16.
The convergence driver and viewers are in the same directory. The optional STEP
file viewer uses MATLAB's PDE Toolbox; the mesher and geometry-package viewer
do not require it. Driver exports go under `outputs/step`.

## Native interface

```sh
./bin/step_mesher_cli mesh \
  --step step_mesher/examples/step_files/test_curved_uv.step \
  --order 2 --edge-order 8 --mesh-fraction 0.10 --refinement 0 \
  --profile rigid --out outputs/step/cli-example \
  --occ-info-exe "$PWD/matlab/+surfsmooth3d/+stepmesher/private/step_mesher_occ_info"
```

The CLI retains equispaced diagnostic samples and CSV output. It does not export
RV geometry; use the MATLAB workflow for GO3/CAD package generation. The native
API also accepts caller-supplied nodes. Omitting or passing a nonexistent helper
path preserves the original Gmsh numerical path. The matching bundled helper
selects original OCC calculations directly in process, without spawning it.
Arbitrary third-party helper executables remain unsupported. MATLAB's unchanged
option defaults select the bundled helper when available.

Original face labels are source group identifiers, not live Gmsh handles after
fragmentation. Current-face ancestry drives Gmsh operations; OCC calculations
retain the original support domains. Metadata reports `native_occ_topology_v1`
and the chosen numerical path. Mesher calls are serialized and clear Gmsh
models, preserving the original lifecycle.

## Projection performance

RV-node projection reuses an initialized OCC projector across consecutive
points on the same bounded B-spline or Bezier face. Complex spline surfaces can
have expensive surface sampling grids; keeping that grid avoids rebuilding it
for every point. Each point still uses OCC's original global projection
algorithm, bounds, tolerance, extrema selection, and output order. Mesh sizing,
quadrature order, and refinement settings are unchanged.

Analytic surfaces and single-point runs use the original projector path. Only
one reusable projector exists at a time, owned by the current face run inside
one `project_nodes` call. A failed projection discards it and retries the
original path, which is then used for the rest of that run. No geometry cache
survives the call or is shared between threads.

`test_projection_reuse` compares XYZ, UV, distances, and fallback counts with
the original fresh-projector calculation, including seams, repeated face
visits, singleton batches, and failed-projection recovery. CTest runs the
Bezier checks. To check another fixture explicitly:

```sh
build/step/test_projection_reuse \
  step_mesher/examples/step_files/wobbly_hairy_torus_10_v2_flat_AP214.step \
  --batch-sizes=1,2,16,64,256 --points-per-face=12 --repetitions=1
```

The MATLAB differential drivers are `tests/step/benchmark_projection_reuse.m`
and `tests/step/compare_projection_reuse.m`. Run baseline and candidate builds
in separate MATLAB processes; the driver verifies the selected package and
loaded MEX. Captures include complete mesh fields, `srcvals`, and area/flux
checks, with strict numerical and byte comparisons across repetitions.

For exhaustive validation of expensive geometries,
`tests/step/instrument_projection_benchmark.py` instruments an archived source
copy under `build/` to capture actual RV projection requests and stage timings.
It refuses to modify the production source. `replay_projection_batches.py`
compares every requested projection using separate baseline and candidate
executables; `compare_projection_replay.m` also compares the uninterrupted
MATLAB mesh with those baseline results. Parallel replay times are diagnostics,
not serial end-to-end speedup measurements.

Recorded validation is in `provenance/projection-reuse-validation.json`: all
270,630 current torus nodes match the original projection results, with a
separate byte comparison of the uninterrupted MATLAB mesh. The report also
records the other fixtures, repeated timings, memory measurements, and their
measurement limits.

## Validation and known limitations

CTest covers the original smoke, native contract and OCC-context tests plus
native/file-import equivalence for every included STEP fixture. It also replays
six frozen CAD operation types and checks five retained-context operations and
CAD measurements after actual helper fragmentation. Portable requests, expected
responses and their source hashes are in `tests/step/fixtures/occ_operations`.
The new helper, typed frozen evaluator and recorded baseline agree byte for byte.
MATLAB tests
are in `tests/step`. The full-working-tree import map is
`provenance/step-imports.json`; extraction results are saved separately under
`build/validation/step`. `docs/V3_NATIVE_IDENTITY_MIGRATION.md` is historical V3
evidence, not a claim that its original paths form the new installation.

This extraction does not add V2 repairs, automatic refinement, automatic order
selection, or new numerical tolerances. The existing `same_support` heuristic,
original-edge searches, inactive reader options, and six-decimal edge-tolerance
conversion are retained. Successful construction, area agreement, and
self-convergence do not certify polynomial-to-CAD accuracy or watertightness.
The original torus's irregular coarse-mesh convergence remains a known
limitation; its regression tests are not an absolute accuracy guarantee.

# SurfSmooth3D

SurfSmooth3D builds high-order triangular surface representations from CAD
geometry and existing scaffolds. It contains the tested surface smoother and
STEP mesher V3, extracted into an independent package without changing their
numerical algorithms. FMM3D is an external dependency; fmm3dbie is not required.

The package supports one-step and two-step smoothing, fixed CAD skeletons,
sigma evaluation and gradients, edge-preserving blending, adaptive smoothing,
adaptive blending, MATLAB GUIs and batch exports. The output `.go3` files can
be read by the existing fmm3dbie package.

## Build

The validated local configuration is Apple Silicon, GNU Fortran 14 and MATLAB
R2026a. A Linux profile is provided but must be validated on Linux before
claiming support there. Native builds do not require MATLAB or the STEP mesher.
MATLAB workflows use the JVM for source hashing; GUI workflows also need
MATLAB graphics. STEP builds require CMake 3.20 or newer and a C++17 compiler.

Install FMM3D separately and select its installation explicitly:

```sh
make lib FMM_PREFIX=/path/to/FMM3D/install
make matlab FMM_PREFIX=/path/to/FMM3D/install
OMP_NUM_THREADS=2 make test-native FMM_PREFIX=/path/to/FMM3D/install
```

`FMM_PREFIX` may contain `lib/libfmm3d.a`, `lib-static/libfmm3d.a`, or
`libfmm3d.a`. The extraction reference is FMM3D revision
`3785618eda092a2234def359227a30e68582080c`. Standalone revision
`d2b5e6e983e1ec915d15330ce74df44221a9ac42` also passed the native and MATLAB
equivalence checks. The default prefix under `build/deps/fmm3d-reference/install`
is a local validation installation, not vendored FMM source.

Override compilers and installation paths with command-line variables or an
untracked `make.inc`. Native and MATLAB objects are kept separate because
their BLAS interfaces differ. The Apple Silicon profile preserves the tested
MATLAB OpenMP and runtime-pinning configuration. Restart MATLAB before
rebuilding or replacing loaded MEX files.

## MATLAB

From the package root:

```matlab
addpath(fullfile(pwd,'matlab'));
root = setup_surfsmooth3d();

opts = struct('fcad',fullfile(root,'tests','fixtures','cylinder', ...
    'cylinder_skeleton.txt'),'nquad',12,'rlam',5,'adapt_sigma',1, ...
    'nrefine',0,'two_stage_smoother',true);
S = surfsmooth3d.multiscale_mesher( ...
    fullfile(root,'tests','fixtures','cylinder','cylinder.gidmsh'),8,opts);
plot(S{1});
```

Use `surfsmooth3d.multiscale_mesher_adaptive` for full adaptive smoothing,
`surfsmooth3d.multiscale_mesher_sigma_eval` for sigma and its gradient, and
`surfsmooth3d.adaptiveblend.solve` for adaptive edge blending. Existing options,
defaults and return conventions are preserved. The one-step path remains the
default of `multiscale_mesher`.

Guarded two-stage and adaptive solves automatically defer Newton targets whose
proposed step is nonfinite or leaves the fixed 10R sphere. After the other
targets converge, recovery scans the original pseudonormal near its initial
position and uses safeguarded Brent iterations to select the closest detected
outward crossing of Phi=1/2. Sources and the sigma field stay fixed. Unresolved
roots or invalid recovered patches stop the solve without returning a mesh.
Set `opts.newton_recovery=false` to restore the previous guarded abort behavior;
the default one-stage path is unaffected. The uniform wrapper also supports
`[S,info] = surfsmooth3d.multiscale_mesher(...)`. Its `info.recovery`, and the
corresponding adaptive/blend metadata, contain counts, per-target stage records,
and a preserved diagnostic path when recovery occurred.

Recovery searches at most eight local sigma on each side with 128 subdivisions
per side. It cannot certify all intersections or restore a feature removed by
smoothing. The separate `upward-tree-reuse` FMM optimization is deferred.

The [recovery validation record](provenance/newton-recovery-validation.json)
contains the regression results, input hashes and controlled timings. Complete
cylinder and intersecting-sphere comparisons retain bitwise output equality
and identical FMM calls, with median overhead below 5%. The hairy-torus check
covers 128 targets on the full fixed CAD sources; full-surface acceptance has
not yet been established. Restart MATLAB after rebuilding the pinned gateways.

The main interactive/batch drivers are in `examples/matlab`:

- `driver_two_stage_edge_smooth_batch.m`
- `driver_two_stage_adaptive_smooth.m`
- `driver_two_stage_adaptive_blend.m`

They use the namespaced class `surfsmooth3d.surfer`, so fmm3dbie's global
`surfer` can remain on the MATLAB path. Exchange `.go3` files or raw arrays
between the packages; old serialized `surfer` objects do not automatically
become the new class.

## Optional STEP mesher

The STEP component requires Gmsh 4.12.2 and OCCT 7.9.3 built together in one
prefix, with the matching native ABI. Its dependency lock and original build
recipe are in `step_mesher`. The optional bootstrap uses Python; the mesher
does not use Python at runtime.

```sh
python3 step_mesher/tools/build_native_dependencies.py \
  --archives-dir /path/to/verified/archives
make step STEP_DEPENDENCY_PREFIX="$PWD/build/deps/cad/install"
make install-step PREFIX="$PWD"
```

The bootstrap is offline by default; `--download` explicitly enables downloads.
Alternatively, supply a prebuilt compatible prefix without running Python.
Set `STEP_MEX=OFF` for a native-only STEP build. See `step_mesher/README.md`
for the retained dependency and runtime checks.

```matlab
mesh = surfsmooth3d.stepmesher.mesh_step('model.step', ...
    struct('order',8,'edge_gll_order',8));
pkg = surfsmooth3d.stepmesher.export_geom_package(mesh,'my_cad_package');
```

The ordinary STEP driver retains p8/GLL8; the lower-level defaults remain
p8/GLL16. No meshing settings have been unified or retuned in this extraction.

## Tests and provenance

```matlab
addpath(fullfile(root,'tests','matlab'));
run_surfsmooth3d_tests('all');
```

Run native tests first to generate the Newton-failure diagnostic fixtures.
Tests write generated results under build/output directories or temporary
folders. The optional large geometry GUI diagnostics accept explicit input
packages; the ordinary regression fixtures are included.

Install native and MATLAB products under a selected prefix with `make install
install-matlab PREFIX=/path/to/install`. Add `make install-step
PREFIX=/path/to/install` when STEP support has been built. STEP application
relocation is tested with its CAD dependency prefix fixed; relocating that
dependency prefix requires rebuilding against the new path.

Read `docs/architecture.md` for component boundaries and file formats,
`docs/validation.md` for measured results, and `provenance` for source hashes,
working-tree changes and import mappings. Numerical improvements should be
made separately from this extraction.

`python3 provenance/verify_sources.py` verifies the recorded extraction hashes
without access to the original repositories. Historical absolute STEP paths
inside frozen fixture metadata are descriptive provenance and are not opened
at runtime.

The original fmm3dbie and STEP development checkouts remain unchanged. This
repository has no configured remote and has not been pushed.

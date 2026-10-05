> Historical upstream V3 evidence, copied unchanged below. Its absolute build paths
> and previous installation results describe the original V3 development environment.
> SurfSmooth3D extraction verification is recorded separately in `build/validation/step`.

# V3 native identity migration

This report records the retained-topology implementation, protected numerical
method, and release evidence. Final installation acceptance is recorded in
`build/native-identity-20260923/STATUS.json` and its referenced publication
record, including the fresh execution of the actual installed driver.

## Scope

One retained OCC STEP import replaces the repeated subprocess imports and
approximate face association. Gmsh imports that same retained topology and
returns the existing tag when an individual retained face is imported again.
The context keeps canonical topology-plus-location identity, first-seen face
ordering, separate oriented occurrences, original edges, and adjacency.

`occ_context.cpp` contains the original helper calculations as typed operations.
`occ_info.cpp` is the standalone serialized adapter. Descriptor-cost face
matching exists only in that compatibility adapter and comparison tests, not
in the production mesher. Bounding boxes remain in the original sizing and
measurement calculations. Original-edge distance searches and `same_support`
remain unchanged.

The mesher keeps source support labels distinct from current Gmsh face tags.
Fragmentation ancestry routes Gmsh evaluations through current live children;
OCC calculations retain the original support domains. Missing/conflicting
identity is an explicit error. Existing numerical fallback behavior is retained.

The six extracted operations retain the old helper-failure boundary: an OCC
operation exception selects the original fallback, as a nonzero helper exit
did previously. A typed native-contract exception keeps identity/runtime errors
fatal. Focused regression tests distinguish successful return, ordinary and
non-standard numerical exceptions, and a native identity error. This boundary
correction does not alter the numerical operation bodies.

The native core serializes mesher calls and retains CAD ownership until Gmsh
cleanup finishes. It verifies the loaded Gmsh and OCC runtime before native
pointer import. Process-global reader settings are initialized to the old
helper's effective defaults and restored afterward. Previously inactive reader
option fields remain inactive, as agreed.

## Protected evidence

The starting V3 sources, complete MATLAB package, native binaries, copied CMake
cache, dependency linkage, driver settings and source/STEP hashes are under
`build/native-identity-20260923/baseline`. Existing user output was cloned there
without replacing the originals. The copied CMake cache refers to an older
source location; all migration builds use new directories.

The independent dependency prefix is `build/native-deps/install`. Exact source
archive hashes are in `dependencies.lock.json`; configure options, compiler and
loaded-library evidence are saved with migration results. Dependencies are
Gmsh 4.12.2 and OCCT 7.9.3, AppleClang 16, arm64, Release, C++17. No Gmsh/OCC
library is loaded from another Astra folder or the Python framework.

## Controlled comparisons

The evidence separates four implementations:

1. Archived original V3 binaries and their original runtime.
2. Unchanged original sources rebuilt with the independent V3 runtime.
3. Extracted helper operations on frozen requests and explicit associations.
4. Complete native identity integration.

Using the same runtime, native integration matches the unchanged implementation
exactly on all 72 compared CSV arrays from the twelve CLI examples in both
numerical modes. All twelve MATLAB p2/GLL8 cases have identical areas and
metadata counters. The three p8/GLL8 fixtures (`test_curved_uv`, UV, torus)
have identical mesh arrays, GLL data and `srcvals`: 102 array comparisons.

All six extracted operations match both the archived and common-runtime helpers
byte-for-byte on 54 frozen requests, including dimensions, ordering, UVs and
operation counters. No enlarged numerical tolerance was needed. Native calls
with frozen MATLAB VR coordinates also match MATLAB on 57 array comparisons.

After preserving the helper-level exception fallback, the final rebuilt runtime
was checked again on all 24 CLI cases, twelve MATLAB cases and the three p8
fixtures. All 246 final-to-previous-candidate comparisons are exact: 72 CLI
arrays, 102 p8 mesh/polynomial arrays, and 72 MATLAB diagnostic comparisons.
The long successful-case studies retain their actual earlier binary provenance;
this final comparison binds their unchanged successful numerical path to the
published build. The new failure boundary is covered separately by native tests.

The dependency change does alter some archived scaffolds, including the sphere
examples and piecewise cap. These differences already occur in unchanged code
rebuilt with the common runtime. They must not be represented as changes made
by the numerical refactor. The three-way matrix retains their counts, areas,
counters and numerical arrays separately.

## Import and fragmentation evidence

For every bundled STEP example, retained-native import and the old Gmsh file
import using the same runtime produce identical sampled support positions,
first/second derivatives, UV bounds, face associations, edge adjacency and edge
lengths. Each face was checked on 81 parameter samples.

Synthetic tests cover located instances, reversed uses, shared topology,
coincident descriptors, distinct trims sharing the same support handle, seams,
fragment ancestry, entity-count stability, repeated native destruction and MEX
unload/reload. The cylinder example deletes source tag 1 and creates current
children 4–9. A test-only observer verified 534 Gmsh surface calls against live
face ancestry.

Gmsh's Boolean builder is not universally non-destructive. The retention audit
found added 2D p-curve representations in UV, `cylinder_bezier_v2`, and the
piecewise cap. Their serialized 3D curves, support surfaces and locations remain
unchanged. All 42 post-fragment frozen operation checks and CAD measurements
remain byte-identical. Other nine examples retain identical serialized BReps.
This is measured regression evidence on these fixtures, not a guarantee for
arbitrary untested CAD. No protected-copy geometry or new numerical algorithm
was introduced.

## Numerical and user-visible behavior

- MATLAB keeps `koorn`/fmm3dbie nodes, fitting, differentiation and quadrature.
- The CLI keeps equispaced diagnostic sampling and the original Gmsh default.
- Selecting the bundled helper enables the original OCC path in process.
- Ordinary driver p8/GLL8 and API p8/GLL16 defaults remain unchanged.
- Manual refinement, all presets, user-selected STEP files and driver settings
  are preserved.
- Solver equations, initial guesses, stopping tolerances, fallback rules,
  continuity-helper construction, GLL interpolation and blending are unchanged.
- No repair allowance, automatic refinement, automatic order selection, new
  accuracy gate, or V2 numerical module is present.

At the C++/CLI boundary, an empty or nonexistent helper selection preserves the
Gmsh numerical path. MATLAB's unchanged default-option function fills an empty
selection with the bundled helper; a nonempty nonexistent path still selects
the old fallback. A custom third-party executable is explicitly rejected rather
than executed. Existing source grouping labels are not live Gmsh handles after
fragmentation.

## Performance and memory

Three-run medians use identical frozen MATLAB VR nodes, p8/GLL8, rigid sizing
and refinement zero. Times are construction seconds:

| Geometry | Original runtime | Unchanged code, common runtime | Native V3 |
|---|---:|---:|---:|
| `test_curved_uv` | 0.504688 | 0.487146 | 0.064256 |
| `UV_curved_bezier` | 3.975811 | 3.407562 | 1.839076 |
| `torus_ellipse` | 65.189135 | 53.090322 | 52.233781 |

The original-runtime timing executable was reconstructed with the original
runtime and compiler optimization settings, with numerical parity checked
against the archived binaries. The common-runtime column isolates the refactor.
No case exceeded the agreed 20% construction-regression investigation threshold.
Convergence/profile jobs were running concurrently, so these are observational
medians under background load, not unloaded benchmark measurements. The small
cases benefit most from removing subprocess startup and repeated imports.

Median process peak resident sizes measured by `/usr/bin/time -l` are:

| Geometry | Original runtime (MiB) | Common runtime (MiB) | Native V3 (MiB) |
|---|---:|---:|---:|
| `test_curved_uv` | 43.89 | 36.78 | 38.75 |
| `UV_curved_bezier` | 47.81 | 41.02 | 45.53 |
| `torus_ellipse` | 60.56 | 62.11 | 59.22 |

These process measurements are not simultaneous aggregate peaks of the parent
and the original helper subprocesses. Raw timings and memory measurements are
retained in `results/performance`, separately from MATLAB processing/export
timings in the geometry-run records.

## Refinement comparisons

All 72 requested runs completed successfully: original and V3, three models,
orders 4/6/8, GLL8, and refinement levels 0–3. Sizing is rigid with mesh fraction
0.1 for the curved and UV examples and 0.05 for the torus. All operation counters
agree across the 36 original/V3 pairs. Full-precision total areas, face areas,
fluxes, scaffold statistics, timings and case-specific binary provenance are
retained in the published evidence under `results/convergence-details` and
the two `*-convergence-merged` directories.

The two non-reference cases eligible for the agreed small-self-error gate pass:

| Case | Original self-error | V3 self-error | Allowed V3 error |
|---|---:|---:|---:|
| `test_curved_uv`, p8/r2 | 3.850176723014788e-13 | 3.848168284452965e-13 | 7.700353446029576e-13 |
| `UV_curved_bezier`, p8/r2 | 1.8703558363547433e-13 | 1.867125515393854e-13 | 3.7407116727094867e-13 |

The complete 36-pair regression report passes. Its other cases record successful
construction and comparative diagnostics; the threshold applies only where the
original self-error is at most 1e-12. The three reference-self comparisons are
zero by definition and are not accuracy tests.

The p8/r3 reference areas are:

| Model | Original | V3 |
|---|---:|---:|
| `test_curved_uv` | 2264.1835280920518 | 2264.1835280920518 |
| `UV_curved_bezier` | 175.96833116287013 | 175.96833116287016 |
| `torus_ellipse` | 10311.872827951158 | 10311.872827951143 |

![Original and V3 refinement comparisons](../build/native-identity-20260923/results/convergence-figures/convergence-comparison.png)

The rigid torus study is irregular in both implementations and does not
demonstrate high-order convergence at these settings. V3 reproduces that
behavior. This migration establishes preservation of the original numerical
method, not absolute torus accuracy or a repair of that method's separate
limitations. No independent CAD-accuracy claim follows from these comparisons.

Long torus cases were run in independent processes in parallel. Complete case
records are assembled with individual source paths and hashes; the original
sequential test jobs were stopped only after their requested p4/p6 results were
saved, to avoid repeating the separately running p8 cases. Earlier records and
logs are preserved.

Both requested torus curvature profiles also construct successfully:

| Profile | Original patches / area | Native V3 patches / area |
|---|---|---|
| `curvature_balanced` | 14152 / 10286.908049768908 | 14152 / 10286.908049768912 |
| `curvature_balanced_local` | 11052 / 10286.908029579465 | 11058 / 10286.908032209798 |

The local-profile difference was isolated with an additional unchanged-code
common-runtime run. It produces 11058 patches and area 10286.908032209798,
exactly matching V3, together with identical flux, face areas and operation
counters. Thus this small scaffold change belongs to the dependency rebuild,
not the identity refactor. The original/V3 relative area difference is
2.5569716853008394e-10. No numerical tolerance or meshing option was changed to
obtain this result.

## Release gates

The native smoke, identity/contract and CAD-context tests pass in the sanitizer,
clean Release/MEX, and final install-prefix builds (3/3 in each). ASan/UBSan
instruments the application core and tests; external Gmsh/OCCT libraries and
MATLAB are not sanitizer-instrumented, and macOS leak detection is disabled.

Clean installation and application relocation pass with the original install
both present and temporarily hidden. Fresh CLI and MATLAB processes loaded the
core and MEX from the intended installation; all 29 loaded Gmsh/OCCT images
resolved to V3's independent dependency prefix. The relocation test moves the
application package while keeping that dependency prefix fixed. It does not
claim that moving the entire dependency tree is supported.

The final target-prefix build is separately installed into a staging root and
tested before publication. Its native core SHA-256 is
`e19867e663daba6ddaa062f78f7a1bb06706dd44de51eb142cc50ea327020ec2`;
its MEX SHA-256 is
`e922f89d08964f25c8d5fda295e49875070ed09d2e7564917e499b12b8a2171c`.

The unchanged ordinary driver completed a staged p8/GLL8 torus run at mesh
fraction 0.05 and refinement zero: 604 patches, 27180 nodes, `surfer`
construction, XYZW/package export, reader round trip, and all six figures.
Following checked publication, the actual V3 driver is run again in a fresh
MATLAB process. Its original output directory is preserved, verification output
is kept separately, and actual core/MEX/dependency paths are checked. The
publication record is marked `verified` only after that final run passes.

The publication preserves rollback copies of the previous complete MATLAB
package, runtime and modified source files. The original fast folder, V2,
earlier baselines, user driver selections and preexisting outputs are retained.

Linux verification remains unperformed: no Linux runner was available in this
macOS session.

The acceptance criterion for old self-errors at most 1e-12 is unchanged:
`E_candidate <= max(2*E_old, 100*eps)`. Self-comparison of the p8/refinement-3
reference is not an accuracy test. CAD reference area and old meshes are not
independent exact references.

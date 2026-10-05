# Extraction validation — 5 October 2026

Validated on macOS 14.7.1 Apple Silicon, GNU Fortran 14.2 and MATLAB R2026a.
FMM3D, BLAS and OpenMP configurations are recorded in `provenance`. Most
numerical regression processes used two OpenMP threads.

## Source preservation

All 14 Fortran smoother files remain byte-identical to the tested fmm3dbie
source. Native STEP sources and headers retain their original computational
content; the only changed C++ source text is the owned-helper installation
path, extended for the MATLAB namespace.

Reversing the mechanical MATLAB namespace substitutions leaves only help
examples, package input/output locations and a diagnostic temporary filename
changed in production MATLAB code. Existing formulas, defaults, tolerances,
refinement decisions and error policies were not edited.

`provenance/source-snapshot.json` captures the original working trees,
including untracked V3 source. Import mappings, a readable mechanical diff and
the final extraction hashes are supplied. `verify_sources.py --source-root
/path/to/original/workspace` additionally checks that the original files remain
unchanged. The generated MWrap C gateway and two MATLAB wrappers regenerate
byte-for-byte with MWrap v1.2.

## Native smoother

Six imported suites passed in each of three builds: optimized reference FMM,
full-core bounds-checked debug, and newer standalone FMM. They cover the
small-distance gradient regression, sigma modes and gradients, Newton guards,
adaptive smoothing, blend sessions, fixed CAD/level-set/sigma sources and
repeated calls.

- Legacy order-eight cylinder outputs at refinement levels zero and one are
  byte-identical to the preserved old-package results.
- Eight geometry outputs using newer FMM3D are byte-identical to reference.
- C-string/Fortran outputs and static/shared-library outputs are identical.
- The full debug build differs from optimized output by at most
  approximately 2.45e-14 in the compared arrays.
- Fixed-source and A-B-A comparisons have exactly zero measured difference.
- An external consumer compiles against the installed native modules/library.

See `provenance/native-validation.json` for individual checks and logs.

## MATLAB and direct extraction equivalence

Separate original-package and extracted-package processes used the same CAD
fixture and settings. Nine result groups are exactly equal: source geometry,
edge data, one-step output, two-stage output, all sigma modes/gradients,
uniform blending, adaptive smoothing, adaptive blending and repeated calls.
The comparison includes coefficients, quadrature weights, refinement maps,
indicators, counts, histories and stopping reasons. All five exported GO3
files are byte-identical.

The same comparison passed again using standalone FMM3D revision d2b5e6e.
Reference-linked MEX binaries were restored afterward.

The pure MATLAB, native-interface, GUI/accuracy and diagnostic groups passed:
18 reloadable sphere batch exports, eight real-CAD batch exports, sigma GUI,
adaptive GUI, adaptive-blend GUI and their existing independent accuracy
checks. The adaptive smoother example used 1769 patches instead of 3872
uniform patches at the deepest level, with zero unresolved patches.

Each of the three smoother gateways also ran alone in a fresh MATLAB process,
with only its intended MEX loaded, followed by `clear mex` and successful
normal shutdown. The existing Mac runtime-pinning profile was enabled.

Two optional GUI diagnostics requiring the large external Multiscale CAD
package were not run; their tests remain available with explicit inputs.
Ordinary Newton-failure plots, snapshots and failed batch exports were tested
using included/generated fixtures. See `provenance/matlab-validation.json`.

## STEP V3

Gmsh 4.12.2 and OCCT 7.9.3 were rebuilt from verified archives into a fresh
common dependency prefix. Native STEP and MATLAB MEX binaries were rebuilt
against that prefix. No copied old build cache or old binary was used for the
new implementation.

- 26/26 native CTest cases passed: import/contract/context checks, all twelve
  STEP fixtures, six frozen operations and five retention checks after
  fragmentation.
- Four original-V3 versus extracted-V3 MATLAB cases each passed 92/92 exact
  numerical comparisons, including GO3 export/reload. Models cover a box,
  curved/UV geometry and a torus.
- All six operation outputs are byte-identical to both current old V3 and
  archived responses.
- Retained-context operations and CAD measurements remain equal after
  fragmentation. BREP serialization acquires a p-curve representation, as in
  the existing behavior; byte equality of that serialization is not claimed.
- The STEP MATLAB smoke test and separately installed CLI/MEX passed.
  Runtime tracing found the new core and new CAD dependencies, with no old
  Astra or Python-framework libraries.

See `provenance/step-validation.json`. Existing geometric accuracy and torus
convergence behavior are preserved; these results establish extraction parity.

## Installed pipeline and BIE interoperability

A separate installation produced a new intersection-of-balls CAD package from
STEP, then ran uniform full/partial smoothing, adaptive full smoothing and
adaptive blending. All exports passed existing validity checks; original CAD
input hashes stayed unchanged. Neither original package was on its MATLAB path.
The fresh STEP mesh had 240 patches, adaptive full smoothing produced 534, and
adaptive blending produced 501. The intentionally depth-limited blend run
retained its explicit tolerance-unmet status.

The unchanged fmm3dbie MATLAB and Fortran readers loaded all five reference
output families. MATLAB arrays agree exactly; native source values agree
exactly, with coefficient differences no larger than 2.67e-15 and weight
differences no larger than 5.56e-17. Both MATLAB class namespaces resolve
correctly in either path order.

A further fresh MATLAB process loaded all six original/new smoother gateways
simultaneously. Both path orders, interleaved A-B-A calls and simultaneously
live blend sessions passed. Closing the original session left the new one
intact. The original BIE layer evaluator retained about 5.63e-10 relative
accuracy, and MATLAB exited normally. All six gateway export tables contain
only the two required MATLAB entry symbols. See
`provenance/coexistence-validation.json`.

Two BIE-only processes consumed saved files with SurfSmooth3D absent from the
MATLAB path. A small sphere check achieved relative error 4.27e-7. The complete
fresh STEP → adaptive smoothing → GO3 → unchanged BIE solver chain achieved
relative error 6.17e-7 and GMRES residual 3.75e-8 or smaller.

The application installation may move while the selected CAD dependency
prefix remains fixed. Moving the latter requires rebuilding against its new
path. Linux, other MATLAB versions, and native co-linking with fmm3dbie in a
single executable are not claimed as validated.

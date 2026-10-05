# Architecture and preserved contracts

## Geometry pipeline

STEP V3 produces a CAD package containing `surface.go3`, `scaffold.gidmsh`,
`edges_gll.txt` and `metadata.txt`. Existing packages or supported scaffold
inputs may be used without building the STEP native component.

The Fortran smoother retains all 14 original source files in
`src/multiscale_mesher`. Its FMM adapter uses tree-building and direct-kernel
interfaces as well as the FMM evaluator, so FMM compatibility is tied to a
tested revision, not simply to the existence of a public `lfmm3d` symbol.

The two-stage algorithm first changes the launch scaffold. Original CAD
quadrature positions, normals, weights and the original sigma construction
remain fixed within the run. Refinement changes output resolution without
redefining that level set. Choosing different sigma settings between runs can
change the surface. A partially blended surface is generally not itself the
level set phi=1/2.

Full adaptive smoothing is implemented in Fortran. Adaptive blending keeps
its existing MATLAB adaptation/blending algorithm and its Fortran session
backend. Extraction does not port either implementation to another language.
Existing state ownership and non-reentrant limitations remain; no new
concurrent-call guarantee is introduced.

## Mathematical utilities

Koornwinder routines and RV tables are copied into `src/special_functions`;
MATLAB basis routines are under `surfsmooth3d.internal`. The coherent original
geometry class and polynomial utility groups are retained. These copies can
later be replaced by a shared special-functions dependency in a separate,
tested change.

Native Legendre/printing routines currently remain provided by external
FMM3D, as in the reference build. Any source copies retained for provenance
are excluded from the compilation list, preventing a second symbol provider.
BLAS/LAPACK adapters and sorting routines are ordinary support utilities.

Native and MATLAB builds use separate object and Fortran module directories.
The native adapter handles the native BLAS integer ABI; MATLAB uses its
64-bit adapter. Numerical integer kinds and gateway signatures are unchanged.

The three MATLAB gateways hide linked implementation symbols. Public MATLAB
names are in `surfsmooth3d`, and MEX filenames have unique package prefixes.
The standalone native libraries retain original Fortran symbols, including
some generic helpers inside the smoother core. Co-linking them with fmm3dbie
in one executable is not supported in this first extraction; shared
mathematical helper ownership must be settled before offering that use case.
Separate native applications communicate through GO3. The isolated MATLAB
gateways are tested together in the same MATLAB process.

The optional `surfacemesh_to_surfer` adapter keeps its original external
Surfacefun/Chebfun requirements. Those packages are not required by ordinary
STEP, smoothing or GO3 workflows.

## GO3 interchange

The text header contains polynomial degree p followed by the number of
triangular patches M. There are M*(p+1)*(p+2)/2 RV nodes, using the reference
triangle with vertices (0,0), (1,0), (0,1). Twelve component-major data blocks
contain position, two parametric tangent vectors and the unit normal.

The degree in the header is p, not p-1. Adaptive h-refinement still uses one
common polynomial degree. Patch order and local-coordinate derivative
conventions are preserved.

GO3 does not contain CAD identities, edge connectivity, units, material
labels, parent maps or convergence certification. Existing sidecars retain
workflow metadata and diagnostics. Reaching maximum depth or a point limit
does not mean the requested adaptive tolerance was achieved.

CAD source packages retain their current parent/child mapping conventions.
The low-level point skeleton has seven columns: position (3), normal (3) and
quadrature weight (1). GLL edge V2 records retain CAD edge IDs, face adjacency,
panel coordinates, unit tangents and arclength weights; the V1 reader remains.

## BIE boundary

SurfSmooth3D exports geometry; fmm3dbie reads it and retains the basis,
quadrature, interpolation, near-search and geometric operations needed by
its solvers. No BIE geometry code was deleted or altered during extraction.
Some formulations, such as RWG-based solvers, also need connectivity that
GO3 alone does not supply. This extraction makes no change to those APIs.

## MATLAB runtime compatibility

The Apple Silicon build retains the existing compatibility adapters: first-use
MEX locking, and runtime retention for the adaptive-blend gateway, while
leaving native session cleanup intact. MATLAB links its own OpenMP runtime;
native examples use the native runtime. This is a platform build profile,
not a numerical change or a claim that every installation needs the workaround.

Restart MATLAB before replacing MEX files. The tests check repeated calls,
session cleanup, clear requests and normal process exit with the profile
enabled. Removal of that profile requires separate validation.

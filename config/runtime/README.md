# Local MATLAB runtime compatibility

The Apple Silicon profile enables the existing MATLAB R2026a / GNU Fortran
workaround through `MATLAB_RUNTIME_PIN=1`. Each adapter includes the ordinary
gateway without changing its numerical calls and locks that MEX after first
use. The adaptive-blend adapter additionally retains the actual linked
libgfortran image with `dladdr` / `dlopen(RTLD_NODELETE)`, preserving its existing
session cleanup and session locks.

This is a platform compatibility option, not a numerical Fortran change or a
universal MATLAB requirement. Other profiles default to
`MATLAB_RUNTIME_PIN=0`. Start a fresh MATLAB process before replacing binaries
or switching builds. For the affected local configuration, changing the option
requires deleting/rebuilding the three MEX files before testing in a fresh
process.

Native builds use their configured OpenMP runtime. The Apple Silicon MATLAB
profile links the external **static** FMM3D archive and MATLAB's libomp; it does
not load an external FMM3D dylib that brings in libgomp. Native and MATLAB
BLAS wrappers and module/object directories are separate.

The imported adapters originate in `tools/fmm3dbie/local/`, recorded in
`provenance/fortran-imports.json`. Only include paths and private identifiers
changed in this extraction.

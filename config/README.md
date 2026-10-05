# Build profiles

The root makefile builds independent native and MATLAB object trees. The native
tree uses `src/support/lapack_wrap.f90` and ordinary 32-bit-integer BLAS/LAPACK.
The MATLAB tree uses `lapack_wrap_64.f90` and MATLAB's 64-bit-integer libraries.
Do not mix these archives or their object files.

`PROFILE=apple-silicon` reproduces the local macOS configuration used during
extraction: GNU Fortran, Accelerate for native BLAS, GNU OpenMP for native
programs, and MATLAB's `libomp` for MEX gateways. Its local shutdown compatibility
option is described in [runtime/README.md](runtime/README.md).
`PROFILE=linux` provides GNU compiler and BLAS/LAPACK defaults; it still requires
validation on the user's Linux toolchain.

Override individual settings on the make command line or in an untracked
`make.inc`. Typical settings are `FC`, `FFLAGS`, `MATLAB_ROOT`, `FMM_PREFIX`,
`NATIVE_BLAS_LIBS`, `OMP`, and `MATLAB_RUNTIME_PIN`. Use a separate `BUILD`
directory for a different compiler or an independent debug build. Effective
compiler flags and link settings are recorded in that directory. Selecting a
different external FMM archive or changing its checksum forces relinking.

The extraction reference is FMM3D revision `3785618`, installed independently at
`build/deps/fmm3d-reference/install`. For another installation, set `FMM_PREFIX`
(expected archive: `lib/libfmm3d.a`); `FMM_STATIC` can specify an archive directly.
The smoother uses FMM tree/list helpers in addition to the public potential
interface, so a replacement FMM version needs the equivalence tests, not merely
a successful link. The FMM archive must contain position-independent code for
shared libraries and MEX builds.

`libsurfsmooth3d.a` contains smoother and copied support objects only; it does not
contain FMM archive members. Native executables link the external FMM archive
explicitly. The native shared library and MATLAB gateways also link that archive
explicitly. External FMM and copied support symbols are private in the native
shared library and MEX gateways. The native shared library preserves all
imported smoother-core Fortran globals, including the core's own `cisurf`
helpers. This preserves the existing native ABI; it does not establish safe
native co-loading with an old full fmm3dbie library that exports the same core.
Do not combine this static archive with another static copy of the same
unprefixed smoother/support routines.

The copies of `legeexps.f` and `prini.f` under `src/special_functions` record their
source provenance but are deliberately excluded from compilation. The selected
external FMM archive supplies those routines, just as in the reference build.
This avoids silently selecting among duplicate providers by archive order.

Run `make test-native` for the imported numerical suites and static/shared
driver equivalence, and `make test-c` for the preserved C-string interface. A
full bounds-checked build can be run independently with, for example:

```sh
OMP_NUM_THREADS=2 make BUILD="$PWD/build/debug" \
  FFLAGS='-fPIC -O0 -g -std=legacy -fcheck=all -fbacktrace' test-native
```

On the extraction machine add `-mno-outline-atomics` to those debug flags to
match the Apple Silicon compatibility setting. Expected array-temporary
warnings from the imported Fortran are not bounds failures.

The optional STEP component uses its own CMake build. `STEP_DEPENDENCY_PREFIX`
selects the common Gmsh/OCCT installation, and `STEP_MEX=OFF` requests a native
build. `make step` configures and builds; `make install-step PREFIX=/path`
installs the already-built STEP component to a selected destination.

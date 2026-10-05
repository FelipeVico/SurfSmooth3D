# GNU/Linux defaults. Override compilers, BLAS and MATLAB_ROOT in make.inc.
NATIVE_BLAS_LIBS = -lblas -llapack
MATLAB_OMP_LIBS = -lgomp
MEX_PLATFORM_FLAGS = -L$(shell dirname "$(shell $(FC) -print-file-name=libgfortran.so)")

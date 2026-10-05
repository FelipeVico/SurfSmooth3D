# Reproduce the locally tested Apple Silicon compiler/runtime configuration.
FC = gfortran
CC = /usr/bin/clang
FFLAGS = -fPIC -O3 -std=legacy -w -mno-outline-atomics
CFLAGS = -fPIC -O3 -arch arm64 -std=c99
NATIVE_BLAS_LIBS = -framework Accelerate
SHARED_EXT = dylib
SHARED_FLAGS = -dynamiclib -fPIC -Wl,-install_name,@rpath/libsurfsmooth3d.dylib
SHARED_EXPORT_FLAGS = -Wl,-exported_symbols_list,$(NATIVE)/exports.txt
MATLAB_ROOT = /Applications/MATLAB_R2026a.app
MATLAB_ARCH = maca64
MEX_EXT = mexmaca64
MATLAB_OMP_LIBS = -L$(MATLAB_ROOT)/bin/maca64 -lomp
MATLAB_RUNTIME_PIN = 1
GFORTRAN_LIB_DIR = $(shell dirname "$(shell $(FC) -print-file-name=libgfortran.dylib)")
MEX_PLATFORM_FLAGS = -L$(GFORTRAN_LIB_DIR) MACOSX_DEPLOYMENT_TARGET=14.0 'CFLAGS=$$CFLAGS -std=c99 -Wno-implicit-function-declaration' 'LDFLAGS=$$LDFLAGS -Wl,-rpath,$(MATLAB_ROOT)/bin/maca64 -Wl,-rpath,$(GFORTRAN_LIB_DIR)'

/* Local Apple Silicon / MATLAB R2026a compatibility adapter.
 * Keep the gateway and its Fortran runtime loaded through thread cleanup.
 * Restart MATLAB before replacing it or switching checkouts.
 * The generated upstream gateway and numerical sources stay unchanged.
 */
#include <dlfcn.h>
extern void _gfortran_st_write(void *);
#define mexFunction surfsmooth3d_imported_gateway
#include "../../matlab/surfsmooth3d_adaptive_blend_routs.c"
#undef mexFunction

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    static int locked = 0;
    if (!locked) {
        /* This session gateway registers mexAtExit. Keep libgfortran mapped
         * independently of MATLAB's final MEX cleanup/unload ordering. */
        Dl_info runtime;
        if (!dladdr((const void *)&_gfortran_st_write, &runtime) ||
            !dlopen(runtime.dli_fname, RTLD_NOW | RTLD_NODELETE)) {
            mexErrMsgIdAndTxt("surfsmooth3d:localRuntime",
                "Cannot retain the local Fortran runtime for thread cleanup.");
        }
        mexLock();
        locked = 1;
    }
    surfsmooth3d_imported_gateway(nlhs, plhs, nrhs, prhs);
}

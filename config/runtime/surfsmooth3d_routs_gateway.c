/* Local Apple Silicon / MATLAB R2026a compatibility adapter.
 * Keep the gateway and its Fortran runtime loaded through thread cleanup.
 * Restart MATLAB before replacing it or switching checkouts.
 * The generated upstream gateway and numerical sources stay unchanged.
 */
#define mexFunction surfsmooth3d_imported_gateway
#include "../../matlab/surfsmooth3d_routs.c"
#undef mexFunction

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    static int locked = 0;
    if (!locked) {
        mexLock();
        locked = 1;
    }
    surfsmooth3d_imported_gateway(nlhs, plhs, nrhs, prhs);
}

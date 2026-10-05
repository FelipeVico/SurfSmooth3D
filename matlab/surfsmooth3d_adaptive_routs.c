/* Independent gateway: no generated gateway or legacy signature changes. */
#include "mex.h"
#include <stdint.h>
#include <math.h>

typedef void (*progress_fn)(const int64_t *, const int64_t *, const int64_t *,
    const int64_t *, const double *);

extern void multiscale_mesher_adaptive_c(const char *, const char *, const char *,
    const int64_t *, const int64_t *, const int64_t *, const int64_t *,
    const double *, const double *, const int64_t *, const int64_t *, int64_t *, progress_fn);

static void progress(const int64_t *pass, const int64_t *patches,
    const int64_t *unresolved, const int64_t *depth, const double *eta)
{
    mexPrintf("Adaptive pass %lld: %lld patches, %lld unresolved, depth %lld\n",
        (long long)*pass, (long long)*patches, (long long)*unresolved, (long long)*depth);
    mexPrintf("  max [position tail, normal tail, position check, normal check] = %.3g %.3g %.3g %.3g\n",
        eta[0],eta[1],eta[2],eta[3]);
    mexEvalString("drawnow limitrate nocallbacks");
}

static double scalar(const mxArray *a, double lo, double hi, int integer)
{
    double x;
    if (!mxIsDouble(a) || mxIsComplex(a) || mxGetNumberOfElements(a) != 1)
        mexErrMsgIdAndTxt("adaptivesmoother:argument", "Expected a real double scalar.");
    x = mxGetScalar(a);
    if (!mxIsFinite(x) || x < lo || x > hi || (integer && x != floor(x)))
        mexErrMsgIdAndTxt("adaptivesmoother:argument", "Adaptive argument is out of range.");
    return x;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int64_t filetype, nquad, p, mode, depth, budget, ier = 0;
    double rlam, tol;
    char *fname, *cad, *root;
    if (nrhs != 11 || nlhs != 1)
        mexErrMsgIdAndTxt("adaptivesmoother:arguments", "Expected eleven inputs and one output.");
    for (int i = 0; i < 3; ++i)
        if (!mxIsChar(prhs[i]) || mxGetM(prhs[i]) != 1 || mxGetN(prhs[i]) == 0)
            mexErrMsgIdAndTxt("adaptivesmoother:path", "Paths must be nonempty character rows.");
    filetype = (int64_t)scalar(prhs[3], 1, 5, 1);
    nquad = (int64_t)scalar(prhs[4], 1, 20, 1);
    p = (int64_t)scalar(prhs[5], 1, 20, 1);
    mode = (int64_t)scalar(prhs[6], 0, 3, 1);
    rlam = scalar(prhs[7], 0, mxGetInf(), 0);
    tol = scalar(prhs[8], 0, 1, 0);
    if (rlam <= 0 || tol <= 0 || tol >= 1)
        mexErrMsgIdAndTxt("adaptivesmoother:argument", "Require rlam > 0 and 0 < tolerance < 1.");
    depth = (int64_t)scalar(prhs[9], 0, 20, 1);
    budget = (int64_t)scalar(prhs[10], 1, 9007199254740991., 1);
    fname = mxArrayToString(prhs[0]);
    cad = mxArrayToString(prhs[1]);
    root = mxArrayToString(prhs[2]);
    multiscale_mesher_adaptive_c(fname,cad,root,&filetype,&nquad,&p,&mode,
        &rlam,&tol,&depth,&budget,&ier,progress);
    mxFree(fname); mxFree(cad); mxFree(root);
    plhs[0] = mxCreateDoubleScalar((double)ier);
}

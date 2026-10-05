#include "occ_context_probe.hpp"
#include <mex.h>
#include <exception>
#include <string>
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs) {
  if (nrhs != 1 || nlhs > 1 || !mxIsChar(prhs[0])) mexErrMsgIdAndTxt("stepmesher:test:Input", "Supply the STEP path.");
  char *raw = mxArrayToString(prhs[0]);
  if (!raw) mexErrMsgIdAndTxt("stepmesher:test:Input", "Could not decode STEP path.");
  std::string path(raw);
  mxFree(raw);
  std::string failure;
  try { run_occ_context_probe(path); } catch (const std::exception &error) { failure = error.what(); }
  catch (...) { failure = "Unrecognized native exception"; }
  if (!failure.empty()) mexErrMsgIdAndTxt("stepmesher:test:NativeIdentity", "%s", failure.c_str());
  if (nlhs) plhs[0] = mxCreateLogicalScalar(true);
}

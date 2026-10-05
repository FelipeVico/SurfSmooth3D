#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>

#ifndef MWF77_RETURN
#define MWF77_RETURN int
#endif

#if defined(MWF77_CAPS)
#define MWF77_multiscale_mesher MULTISCALE_MESHER_UNIF_REFINE
#elif defined(MWF77_UNDERSCORE1)
#define MWF77_multiscale_mesher multiscale_mesher_unif_refine_ 
#elif defined(MWF77_UNDERSCORE0)
#define MWF77_multiscale_mesher multiscale_mesher_unif_refine 
#else
#define MWF77_multiscale_mesher multiscale_mesher_unif_refine__
#endif

void multiscale_mesher_unif_refine_cfname_(char *, int64_t*, int64_t*, char *, int64_t*, int64_t *, int64_t *, int64_t *, double *, char *, int64_t *);
void f2cstr_(char *);

#ifdef __cplusplus
extern "C"
#endif

int main(int argc, char **argv)
{
  if (argc != 4) {
    fprintf(stderr, "Usage: %s scaffold.gidmsh cad.txt output_root\n", argv[0]);
    return 2;
  }
  char *filenamein;
  filenamein = argv[1];

  char *filenameout;
  filenameout = argv[3];

  int64_t ifiletype = 3;
  int64_t norder_skel = 16;
  int64_t norder_smooth = 8;
  int64_t nrefine = 0;
  int64_t adapt_flag = 1;
  double rlam = 5.0;
  int64_t ier = 0;
  int64_t ifcad = 1;
  char *fcad;
  fcad = argv[2];
  

  multiscale_mesher_unif_refine_cfname_(filenamein, &ifiletype, &ifcad, fcad, &norder_skel, &norder_smooth, 
     &nrefine, &adapt_flag, &rlam, filenameout, &ier);
  if (ier != 0) {
    fprintf(stderr, "Smoother C interface returned error %lld\n", (long long)ier);
    return 1;
  }
  return 0;

}

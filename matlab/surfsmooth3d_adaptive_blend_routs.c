/* Private, explicitly owned sessions; legacy gateways are not changed. */
#include "mex.h"
#include <stdint.h>
#include <string.h>
#include <math.h>

extern void adaptive_blend_open(const char *,const char *,const char *,const int64_t *,
    const int64_t *,const int64_t *,const double *,const int64_t *,int64_t *,int64_t *);
extern void adaptive_blend_close(const int64_t *,int64_t *);
extern void adaptive_blend_close_all(void);
extern void adaptive_blend_size(const int64_t *,int64_t *,int64_t *,double *,int64_t *);
extern void adaptive_blend_get(const int64_t *,const int64_t *,const int64_t *,double *,double *);
extern void adaptive_blend_refine(const int64_t *,const int64_t *,const int64_t *,const int64_t *,int64_t *);
extern void adaptive_blend_project(const int64_t *,const int64_t *,const int64_t *,const double *,double *,int64_t *);
extern void adaptive_blend_sigma(const int64_t *,const int64_t *,const double *,double *,int64_t *);

static int registered=0;
static void cleanup(void) { adaptive_blend_close_all(); }
static void require(int condition, const char *message)
{
    if (!condition) mexErrMsgIdAndTxt("adaptiveblend:gateway", "%s",message);
}
static double scalar(const mxArray *a,double lo,double hi,int integer)
{
    require(mxIsDouble(a) && !mxIsComplex(a) && !mxIsSparse(a) && mxGetNumberOfElements(a)==1,
        "Expected a full real double scalar.");
    double x=mxGetScalar(a);
    require(mxIsFinite(x) && x>=lo && x<=hi && (!integer || floor(x)==x),"Scalar out of range.");
    return x;
}
static void matrix(const mxArray *a,mwSize rows)
{
    require(mxIsDouble(a) && !mxIsComplex(a) && !mxIsSparse(a) && mxGetNumberOfDimensions(a)==2 &&
        mxGetM(a)==rows && mxGetN(a)>0,"Invalid real matrix shape.");
    double *x=mxGetPr(a);
    for (mwSize k=0;k<mxGetNumberOfElements(a);++k) require(mxIsFinite(x[k]),"Nonfinite input.");
}
static int64_t *indices(const mxArray *a,int64_t count,int sorted)
{
    require(mxIsDouble(a) && !mxIsComplex(a) && !mxIsSparse(a) && mxGetNumberOfElements(a)>0,
        "Expected nonempty real indices.");
    double *x=mxGetPr(a);
    mwSize n=mxGetNumberOfElements(a);
    for (mwSize j=0;j<n;++j) {
        require(mxIsFinite(x[j]) && x[j]>=1 && x[j]<=count && floor(x[j])==x[j],"Invalid leaf index.");
        if (sorted && j) require(x[j]>x[j-1],"Refinement IDs must be sorted and unique.");
    }
    int64_t *out=mxMalloc(n*sizeof(int64_t));
    for (mwSize j=0;j<n;++j) out[j]=(int64_t)x[j];
    return out;
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    char command[32];
    require(nrhs>=1 && mxIsChar(prhs[0]) && mxGetM(prhs[0])==1,"Expected a command.");
    require(mxGetString(prhs[0],command,sizeof(command))==0,"Invalid command.");
    if (!registered) { mexAtExit(cleanup); registered=1; }
    int64_t id=0,ier=0,n=0,np=0;
    double sphere[4];
    if (!strcmp(command,"open")) {
        require(nrhs==9 && nlhs==2,"open expects eight arguments and two outputs.");
        for (int j=1;j<=3;++j) require(mxIsChar(prhs[j]) && mxGetM(prhs[j])==1 && mxGetN(prhs[j])>0,"Invalid path.");
        int64_t nq=(int64_t)scalar(prhs[4],1,20,1),p=(int64_t)scalar(prhs[5],1,20,1);
        int64_t mode=(int64_t)scalar(prhs[6],0,3,1),budget=(int64_t)scalar(prhs[8],1,9007199254740991.,1);
        double rlam=scalar(prhs[7],0,mxGetInf(),0);
        require(rlam>0,"rlam must be positive.");
        char *fname=mxArrayToString(prhs[1]),*cad=mxArrayToString(prhs[2]),*root=mxArrayToString(prhs[3]);
        adaptive_blend_open(fname,cad,root,&nq,&p,&mode,&rlam,&budget,&id,&ier);
        mxFree(fname); mxFree(cad); mxFree(root);
        if (!ier) mexLock();
        plhs[0]=mxCreateDoubleScalar((double)id); plhs[1]=mxCreateDoubleScalar((double)ier);
        return;
    }
    require(nrhs>=2,"Missing session handle.");
    id=(int64_t)scalar(prhs[1],1,9007199254740991.,1);
    if (!strcmp(command,"close")) {
        require(nrhs==2 && nlhs==0,"close takes only a handle and no outputs.");
        adaptive_blend_close(&id,&ier);
        if (!ier) mexUnlock();
        return;
    }
    adaptive_blend_size(&id,&n,&np,sphere,&ier);
    require(!ier,"Session is closed or invalid.");
    if (!strcmp(command,"get")) {
        require(nrhs==2 && nlhs==3,"get returns values, leaf maps and enclosing sphere.");
        plhs[0]=mxCreateDoubleMatrix(12,np*n,mxREAL);
        plhs[1]=mxCreateDoubleMatrix(9,n,mxREAL);
        plhs[2]=mxCreateDoubleMatrix(4,1,mxREAL);
        adaptive_blend_get(&id,&n,&np,mxGetPr(plhs[0]),mxGetPr(plhs[1]));
        memcpy(mxGetPr(plhs[2]),sphere,4*sizeof(double));
    } else if (!strcmp(command,"refine")) {
        require(nrhs==4 && nlhs==1,"refine expects IDs and budget, returns status.");
        int64_t budget=(int64_t)scalar(prhs[3],1,9007199254740991.,1);
        int64_t *ids=indices(prhs[2],n,1),nt=(int64_t)mxGetNumberOfElements(prhs[2]);
        adaptive_blend_refine(&id,&nt,ids,&budget,&ier); mxFree(ids);
        plhs[0]=mxCreateDoubleScalar((double)ier);
    } else if (!strcmp(command,"project")) {
        require(nrhs==4 && nlhs==2,"project expects IDs and reference coordinates, returns values and status.");
        matrix(prhs[3],2);
        int64_t nt=(int64_t)mxGetN(prhs[3]);
        require(mxGetNumberOfElements(prhs[2])==(mwSize)nt,"One leaf ID per reference point is required.");
        int64_t *ids=indices(prhs[2],n,0);
        plhs[0]=mxCreateDoubleMatrix(12,nt,mxREAL);
        adaptive_blend_project(&id,&nt,ids,mxGetPr(prhs[3]),mxGetPr(plhs[0]),&ier); mxFree(ids);
        plhs[1]=mxCreateDoubleScalar((double)ier);
    } else if (!strcmp(command,"sigma")) {
        require(nrhs==3 && nlhs==2,"sigma expects targets, returns sigma/gradient and status.");
        matrix(prhs[2],3);
        int64_t nt=(int64_t)mxGetN(prhs[2]);
        plhs[0]=mxCreateDoubleMatrix(4,nt,mxREAL);
        adaptive_blend_sigma(&id,&nt,mxGetPr(prhs[2]),mxGetPr(plhs[0]),&ier);
        plhs[1]=mxCreateDoubleScalar((double)ier);
    } else require(0,"Unknown adaptive blend command.");
}

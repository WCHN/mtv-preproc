#include "mex.h"
#include <math.h>

float samp0(mwSize dm[], float f[], float x, float y, float z)
{
    mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
    float f000, f100, f010, f110, f001, f101, f011, f111;
    float dx1, dx2, dy1, dy2, dz1, dz2;

    ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
    iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
    iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;
    ix1  = ix+1;
    iy1  = iy+1;
    iz1  = iz+1;

    if (iz>=0 && iz<dm[2])
    {
        mwSize tmpz  = dm[1]*iz;
        if (iy>=0 && iy<dm[1])
        {
            mwSize tmpy  = dm[0]*(iy + tmpz);
            f000 = (ix >=0 && ix <dm[0]) ? f[ix +tmpy] : 0.0;
            f100 = (ix1>=0 && ix1<dm[0]) ? f[ix1+tmpy] : 0.0;
        }
        else
            f000 = f100 = 0.0;

        if (iy1>=0 && iy1<dm[1])
        {
            mwSize tmpy  = dm[0]*(iy1 + tmpz);
            f010 = (ix >=0 && ix <dm[0]) ? f[ix +tmpy] : 0.0;
            f110 = (ix1>=0 && ix1<dm[0]) ? f[ix1+tmpy] : 0.0;
        }
        else
            f010 = f110 = 0.0;
    }
    else
        f000 = f100 = f010 = f110 = 0.0;

    if (iz1>=0 && iz1<dm[2])
    {
        mwSize tmpz  = dm[1]*iz1;
        if (iy>=0 && iy<dm[1])
        {
            mwSize tmpy  = dm[0]*(iy + tmpz);
            f001 = (ix >=0 && ix <dm[0]) ? f[ix +tmpy] : 0.0;
            f101 = (ix1>=0 && ix1<dm[0]) ? f[ix1+tmpy] : 0.0;
        }
        else
            f001 = f101 = 0.0;

        if (iy1>=0 && iy1<dm[1])
        {
            mwSize tmpy  = dm[0]*(iy1 + tmpz);
            f011 = (ix >=0 && ix <dm[0]) ? f[ix +tmpy] : 0.0;
            f111 = (ix1>=0 && ix1<dm[0]) ? f[ix1+tmpy] : 0.0;
        }
        else
            f011 = f111 = 0.0;
    }
    else
        f001 = f101 = f011 = f111 = 0.0;

    return( ((f000*dx2 + f100*dx1)*dy2 + (f010*dx2 + f110*dx1)*dy1)*dz2
          + ((f001*dx2 + f101*dx1)*dy2 + (f011*dx2 + f111*dx1)*dy1)*dz1 );
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float *f, *Y, *wf, *x, *y, *z;
    mwSize nd, i, mm;
    mwSize dmf[4];
    const mwSize *dmyp;

    if (nrhs == 0) mexErrMsgTxt("Incorrect usage");
    if (nrhs == 2)
    {
        if (nlhs > 1) mexErrMsgTxt("Only 1 output argument required");
    }
    else
        mexErrMsgTxt("Two input arguments required");

    for(i=0; i<nrhs; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>3) mexErrMsgTxt("Wrong number of dimensions.");
    dmf[0] = dmf[1] = dmf[2] = 1;
    for(i=0; i<nd; i++)
        dmf[i] = mxGetDimensions(prhs[0])[i];

    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dmyp = mxGetDimensions(prhs[1]);
    if (dmyp[3]!=3)
        mexErrMsgTxt("Incompatible dimensions.");

    plhs[0] = mxCreateNumericArray(3,dmyp, mxSINGLE_CLASS, mxREAL);

    f = (float *)mxGetPr(prhs[0]);
    Y = (float *)mxGetPr(prhs[1]);
    wf= (float *)mxGetPr(plhs[0]);

    mm  = dmyp[0]*dmyp[1]*dmyp[2];
    x   = Y;
    y   = Y+mm;
    z   = Y+mm*2;
    for (i=0; i<mm; i++)
    {
        wf[i] = samp0(dmf, f, x[i]-1.0, y[i]-1.0, z[i]-1.0);
    }
}


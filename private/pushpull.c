#include <math.h>
#include "mex.h"
#include "shoot_boundary.h"
#include "fast.h" // fast approximations of exp/log/erf/...

//-------------------------------------------------------------------------
// Main function
//-------------------------------------------------------------------------

/* TODO
 * - Remove out of bound checks
 * - Include other selection profiles (rect, sinc)
 */

/** Pull or push an image according to a deformation and a window
 *
 * @param dm0   [mwSize 3]          Reference dim    (pull:in / push:out)
 * @param  m1   [mwSize]            Nb pulled voxels (pull:out / push:in)
 * @param  n    [mwSize]            Nb features      (4th dimension)
 * @param Psi   [float m1*3]        Deformation
 * @param  F0   [float prod(dm0)*n] Reference volume (pull:in / push:out)
 * @param  S0   [float]             Count volume     (push:out)
 * @param  F1   [float m1*n]        Pulled volume    (pull:out / push:in)
 * @param code  [uint]              (0=push|1=pull|2=pushc|3=pullc)
 * @param Jac   [float m1*3*3]      Jacobian of the deformation
 * @param mJ    [mwSize]            Number of Jacobians (1 or m1)
 * @param func  [*() 3]             Selection window (&gauss/&rec/&sinc)
 */
#define TINY 5e-2f
static void jpushpull(mwSize dm0[], mwSize m1, mwSize n, 
                      float Psi[], float F0[], float S0[], float F1[], 
                      unsigned int code, float Jac[], mwSize mJ)
                      /* (float (*func[])(float,float,float) */
{
    /* Pointers into input/output arrays */
    mwSize  m0 = dm0[0]*dm0[1]*dm0[2],  /* Number of reference voxels */
            oy = dm0[0],                /* Offsets between lines */
            oz = dm0[0]*dm0[1];         /* Offsets between slices */
    float  *px  = Psi;                  /* Deformation in x */
    float  *py  = Psi+m1;               /*                y */
    float  *pz  = Psi+m1*2;             /*                z */
    float  *jxx = Jac;                  /* Jacobian in xx */
    float  *jyx = Jac+mJ;               /*             yx */
    float  *jzx = Jac+mJ*2;             /*             zx */
    float  *jxy = Jac+mJ*3;             /*             xy */
    float  *jyy = Jac+mJ*4;             /*             yy */
    float  *jzy = Jac+mJ*5;             /*             zy */
    float  *jxz = Jac+mJ*6;             /*             xz */
    float  *jyz = Jac+mJ*7;             /*             yz */
    float  *jzz = Jac+mJ*8;             /*             zz */
    float  NaN  = mxGetNaN();
    
    
    /* Stuff needed inside the loop */
    float  Txx, Tyy, Tzz,       /* Inverse covariance */
           Txy, Txz, Tyz,
           Tyx, Tzx, Tzy;
    float  limx, limy, limz;    /* Bounding box in ref space */
    float  norm;                /* Kernel normalisation */
    
    /* Constants (should probably be static outside of the function) */
    const float  isig2  = 8.*log(2.);           /* Gaussian sd so that FWHM = 1 px */
    const float  isig   = sqrt(isig2);
    const float  sig    = 1./isig;
    float  pinorm = 1./sqrt(8*M_PI*M_PI*M_PI);  /* Part of kernel normalisation */
           pinorm = pinorm*pinorm*pinorm;       /* 3D so power 3 */
    
    int voxJ = (mJ != 1);
    if (!voxJ)
    /* Invert Jacobian if same for all voxels */
    {
        /* Jacobian at point i */
        float  Jxx = *(jxx),
               Jyy = *(jyy),
               Jzz = *(jzz),
               Jxy = *(jxy),
               Jxz = *(jxz),
               Jyz = *(jyz),
               Jyx = *(jyx),
               Jzx = *(jzx),
               Jzy = *(jzy);
        
        /* Precompute some values */
        float  JxxJyx = Jxx*Jyx,
               JxxJzx = Jxx*Jzx,
               JyyJxy = Jyy*Jxy,
               JyyJzy = Jyy*Jzy,
               JzzJxz = Jzz*Jxz,
               JzzJyz = Jzz*Jyz,
               JxzJyz = Jxz*Jyz,
               JyxJzx = Jyx*Jzx,
               JxyJzy = Jxy*Jzy;
        
        /* Covariance = J*J'+I */
        float  Cxx = Jxx*Jxx + Jxy*Jxy + Jxz*Jxz + 1,
               Cxy = JxxJyx  + JyyJxy  + JxzJyz,
               Cxz = JxxJzx  + JxyJzy  + JzzJxz,
               Cyx = JxxJyx  + JyyJxy  + JxzJyz,
               Cyy = Jyx*Jyx + Jyy*Jyy + Jyz*Jyz + 1,
               Cyz = JyxJzx  + JyyJzy  + JzzJyz,
               Czx = JxxJzx  + JxyJzy  + JzzJxz,
               Czy = JyxJzx  + JyyJzy  + JzzJyz,
               Czz = Jzx*Jzx + Jzy*Jzy + Jzz*Jzz + 1;
        
        /* Invert covariance */
        float idt = 1.0/(Cxx*Cyy*Czz + 2*Cxy*Cxz*Cyz
                    -Cxx*Cyz*Cyz-Cyy*Cxz*Cxz-Czz*Cxy*Cxy);
        float scl = idt*isig2;
        Txx = scl*(Cyy*Czz-Cyz*Cyz);
        Tyx = scl*(Cxz*Cyz-Cxy*Czz);
        Tzx = scl*(Cxy*Cyz-Cxz*Cyy);
        Txy = scl*(Cxz*Cyz-Cxy*Czz);
        Tyy = scl*(Cxx*Czz-Cxz*Cxz);
        Tzy = scl*(Cxy*Cxz-Cxx*Cyz);
        Txz = scl*(Cxy*Cyz-Cxz*Cyy);
        Tyz = scl*(Cxy*Cxz-Cxx*Cyz);
        Tzz = scl*(Cxx*Cyy-Cxy*Cxy);
        
        /* Normalising constant = (2*pi*sig2)^(-3/2) * det(T)^(1/2) */
        norm = sqrt(idt)*pinorm*isig2*isig;
        
        /* Bounding box of contributing points in ref space */
        limx = limy = limz = 0;
        limx = fmaxf(limx, fabs(Jxx));
        limx = fmaxf(limx, fabs(Jxy));
        limx = fmaxf(limx, fabs(Jxz));
        limy = fmaxf(limy, fabs(Jyx));
        limy = fmaxf(limy, fabs(Jyy));
        limy = fmaxf(limy, fabs(Jyz));
        limz = fmaxf(limz, fabs(Jzx));
        limz = fmaxf(limz, fabs(Jzy));
        limz = fmaxf(limz, fabs(Jzz));
        limx += 1;
        limy += 1;
        limz += 1;
        limx *= 3*sig;
        limy *= 3*sig;
        limz *= 3*sig;
//         mexPrintf("Limits: %f %f %f\n", limx, limy, limz);
    }
    
    /* Precompute exponential */
    float  pexp[1024];
    float  scl = -0.5*9./1023.;
    for(mwSize k=0; k<1024; ++k) pexp[k] = exp(scl*(float)k)*norm;
    scl = 1023./9;
    
    float *pf1 = F1;
    for (mwSize i=0; i<m1; ++i, ++pf1)    /* loop over pulled voxels */
    {
        float  x, y, z;     /* Coordinates in reference volume */
        x = *(px++)-1.0f;   /* Subtract 1 because of MATLAB indexing */
        y = *(py++)-1.0f;
        z = *(pz++)-1.0f;
        
        if (((code & 2)==2 && mxIsFinite(x) && mxIsFinite(y) && mxIsFinite(z))
            || ((x>=-TINY) && (x<=(float)(dm0[0])-1.0f+TINY)
            &&  (y>=-TINY) && (y<=(float)(dm0[1])-1.0f+TINY)
            &&  (z>=-TINY) && (z<=(float)(dm0[2])-1.0f+TINY)))
        /* If (pullc/pushc) OR (pull/push AND inside bounds) */
        {
            if (voxJ)
            /* Invert Jacobian if voxel-specific */
            {
                mexErrMsgTxt("Not implemented yet");
            }
            
            /* Voxels inside the bounding box */
            mwSignedIndex  ixmax = floor(x+limx),
                           iymax = floor(y+limy),
                           izmax = floor(z+limz),
                           ixmin = ceil(x-limx),
                           iymin = ceil(y-limy),
                           izmin = ceil(z-limz);
            mwSize         lx    = ixmax - ixmin + 1,
                           ly    = iymax - iymin + 1,
                           lz    = izmax - izmin + 1;
            
            /* Create lookups of voxel locations - for coping with edges */
            mwSize  oox[128], ooy[128], ooz[128];
            for(mwSize k=0; k<lx; ++k) oox[k] = bound(ixmin+k, dm0[0]);
            for(mwSize k=0; k<ly; ++k) ooy[k] = bound(iymin+k, dm0[1])*oy;
            for(mwSize k=0; k<lz; ++k) ooz[k] = bound(izmin+k, dm0[2])*oz;
            float  ddx[128], ddy[128], ddz[128];
            for(mwSize k=0; k<lx; ++k) ddx[k] = x-(float)(k+ixmin);
            for(mwSize k=0; k<ly; ++k) ddy[k] = y-(float)(k+iymin);
            for(mwSize k=0; k<lz; ++k) ddz[k] = z-(float)(k+izmin);
            
            if ((code&1)==1)
            /* Pull */
            {
                            
                /* Loop over contributing voxels in ref space */
                for (mwSize iz=0; iz<lz;)
                {
                    float dz    = ddz[iz];
                    float *pf0z = F0 + ooz[iz++];
                    for (mwSize iy=0; iy<ly;)
                    {
                        float dy    = ddy[iy];
                        float *pf0y = pf0z + ooy[iy++];
                        for (mwSize ix=0; ix<lx;)
                        {
                            float dx    = ddx[ix];
                            float *pf0x = pf0y + oox[ix++];
                            
//                             float Tdx  = Txx * dx + Txy * dy + Txz * dz;
//                             float Tdy  = Tyx * dx + Tyy * dy + Tyz * dz;
//                             float Tdz  = Tzx * dx + Tzy * dy + Tzz * dz;
//                             float Td2  = dx * Tdx + dy * Tdy + dz * Tdz;
                            float Td2 = dx * (Txx * dx + Txy * dy + Txz * dz)
                                      + dy * (Tyx * dx + Tyy * dy + Tyz * dz)
                                      + dz * (Tzx * dx + Tzy * dy + Tzz * dz);
                            
                            if (Td2 < 9.)
                            {
//                                 float w    = fasterexp(-0.5*Td2)*norm;
                                float w = pexp[(mwSize)floor(Td2*scl)];
                            
                                float *pf1k = pf1;
                                for (mwSize k=0; k<n; ++k, pf0x += m0, pf1k += m1)
                                {
//                                     if (pf0x-F0 <0 || pf0x-F0 >= m0*n)
//                                         mexErrMsgTxt("Out of bound! (pull acc)");
                                    (*pf1k) += w * (*pf0x);
                                }
                            }
                        }
                    }
                }
            }
            else
            /* Push */
            {
                /* Loop over contributing voxels in ref space */
                for (mwSize iz=0; iz<lz;)
                {
                    float dz    = ddz[iz];
                    float *pf0z = F0 + ooz[iz];
                    float *ps0z = S0 + ooz[iz++];
                    for (mwSize iy=0; iy<ly;)
                    {
                        float dy    = ddy[iy];
                        float *pf0y = pf0z + ooy[iy];
                        float *ps0y = ps0z + ooy[iy++];
                        for (mwSize ix=0; ix<lx;)
                        {
                            float dx    = ddx[ix];
                            float *pf0x = pf0y + oox[ix];
                            float *ps0x = ps0y + oox[ix++];
                            
//                             float Tdx  = Txx * dx + Txy * dy + Txz * dz;
//                             float Tdy  = Tyx * dx + Tyy * dy + Tyz * dz;
//                             float Tdz  = Tzx * dx + Tzy * dy + Tzz * dz;
//                             float Td2  = dx * Tdx + dy * Tdy + dz * Tdz;
                            float Td2 = dx * (Txx * dx + Txy * dy + Txz * dz)
                                      + dy * (Tyx * dx + Tyy * dy + Tyz * dz)
                                      + dz * (Tzx * dx + Tzy * dy + Tzz * dz);
                            
                            if (Td2 < 9.)
                            {
//                                 float w    = fasterexp(-0.5*Td2)*norm;
                                float w = pexp[(mwSize)floor(Td2*scl)];

                                float *pf1k = pf1;
                                for (mwSize k=0; k<n; ++k, pf0x += m0, pf1k += m1)
                                {
//                                     if (pf0x-F0 <0 || pf0x-F0 >= m0*n)
//                                         mexErrMsgTxt("Out of bound! (push F0).");
//                                     if (pf1k-F1 <0 || pf1k-F1 >= m1*n)
//                                         mexErrMsgTxt("Out of bound! (push F1).");
                                    (*pf0x) += w * (*pf1k);
                                }
                                /* Count image */
                                if (S0!=0)
                                {
//                                     if (ps0x-S0 <0 || ps0x-S0 >= m0)
//                                         mexErrMsgTxt("Out of bound! (push S0).");
                                    (*ps0x) += w;
                                }
                            }
                        }
                    }
                }
            }
        }
        else if ((code&1)==1)
        /* If (pull/push AND out of bounds) */
        {
            float *pf1k = pf1;
            for(mwSize k=0; k<n; ++k, pf1k += m1)
            {
//                 if (pf1k-F1 <0 || pf1k-F1 >= m1*n)
//                     mexErrMsgTxt("Out of bound! (NaN).");
                (*pf1k) = NaN;
            }
        }
    }
}

//-------------------------------------------------------------------------
// Wrappers
//-------------------------------------------------------------------------

static void jpullc(mwSize dm0[], mwSize m1, mwSize n, float Psi[],  float F0[], /*@out@*/ float F1[], float Jac[], mwSize mJ)
{
    jpushpull(dm0, m1, n, Psi, F0, 0, F1, 3, Jac, mJ);
}

static void  jpull(mwSize dm0[], mwSize m1, mwSize n, float Psi[],  float F0[], /*@out@*/ float F1[], float Jac[], mwSize mJ)
{
    jpushpull(dm0, m1, n, Psi, F0, 0, F1, 1, Jac, mJ);
}

static void jpushc(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F1[], /*@out@*/ float F0[], /*@null@@out@*/ float S0[], float Jac[], mwSize mJ)
{
    jpushpull(dm0, m1, n, Psi, F0, S0, F1, 2, Jac, mJ);
}

static void  jpush(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F1[], /*@out@*/ float F0[], /*@null@@out@*/ float S0[], float Jac[], mwSize mJ)
{
    jpushpull(dm0, m1, n, Psi, F0, S0, F1, 0, Jac, mJ);
}

//-------------------------------------------------------------------------
// Parse MATLAB arguments
//-------------------------------------------------------------------------

static void jpull_mexFunction(int flag, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize nd, i, dm0[4], dm1[4];
    const mwSize *dmpsi, *dmjac;

    if (nrhs == 0) mexErrMsgTxt("Incorrect usage");
    if (nrhs == 3)
    {
        if (nlhs > 1) mexErrMsgTxt("Only 1 output argument required");
    }
    else
        mexErrMsgTxt("Three input arguments required");

    for(i=0; i<nrhs; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");

    /* Dimensions of image to resample */
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    dm0[0] = dm0[1] = dm0[2] = dm0[3] = 1;
    for(i=0; i<nd; i++)
        dm0[i] = mxGetDimensions(prhs[0])[i];

    /* Dimensions of deformation field and resulting output volume */
    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dmpsi = mxGetDimensions(prhs[1]);
    if (dmpsi[3]!=3)
        mexErrMsgTxt("Incompatible dimensions.");
    
    /* Dimensions of Jacobian field */
    nd = mxGetNumberOfDimensions(prhs[2]);
    if (nd!=5) mexErrMsgTxt("Wrong number of dimensions.");
    dmjac = mxGetDimensions(prhs[2]);
    if ((dmjac[0]!=dmpsi[0] || dmjac[1]!=dmpsi[1] || dmjac[2]!=dmpsi[2]) &&
        (dmjac[0]*dmjac[1]*dmjac[2]>1))
        mexErrMsgTxt("Deformation and Jacobian fields should have the same size.");
    if (dmjac[3]!=3 || dmjac[4]!=3)
        mexErrMsgTxt("Incompatible dimensions.");

    /* Dimensions of output volumes */
    dm1[0] = dmpsi[0];
    dm1[1] = dmpsi[1];
    dm1[2] = dmpsi[2];
    dm1[3] = dm0[3];  /* Number of output volumes */
    plhs[0] = mxCreateNumericArray(4,dm1, mxSINGLE_CLASS, mxREAL);
    
    if ((flag&2)==2)
        jpull (dm0, dm1[0]*dm1[1]*dm1[2], dm0[3], (float *)mxGetPr(prhs[1]), (float *)mxGetPr(prhs[0]), (float *)mxGetPr(plhs[0]), (float *)mxGetPr(prhs[2]), dmjac[0]*dmjac[1]*dmjac[2]);
    else
        jpullc(dm0, dm1[0]*dm1[1]*dm1[2], dm0[3], (float *)mxGetPr(prhs[1]), (float *)mxGetPr(prhs[0]), (float *)mxGetPr(plhs[0]), (float *)mxGetPr(prhs[2]), dmjac[0]*dmjac[1]*dmjac[2]);
}

static void jpush_mexFunction(int flag, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float *S0;
    mwSize nd, i, dm0[4], dm1[4];
    const mwSize *dmpsi, *dmjac;

    if ((nrhs != 3) && (nrhs != 4))
        mexErrMsgTxt("Three or four input arguments required");
    if (nlhs  > 2) mexErrMsgTxt("Up to two output arguments required");
    
    for(i=0; i<3; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) ||
              mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");

    /* Dimensions of input volumes */
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    dm1[0] = dm1[1] = dm1[2] = dm1[3] = 1;
    for(i=0; i<nd; i++)
        dm1[i] = mxGetDimensions(prhs[0])[i];

    /* Dimensions of deformation */
    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dmpsi = mxGetDimensions(prhs[1]);
    if (dmpsi[0]!=dm1[0] || dmpsi[1]!=dm1[1] || dmpsi[2]!=dm1[2] || dmpsi[3]!=3)
        mexErrMsgTxt("Incompatible dimensions.");
    
    /* Dimensions of Jacobian field */
    nd = mxGetNumberOfDimensions(prhs[2]);
    if (nd!=5) mexErrMsgTxt("Wrong number of dimensions.");
    dmjac = mxGetDimensions(prhs[2]);
    if ((dmjac[0]!=dmpsi[0] || dmjac[1]!=dmpsi[1] || dmjac[2]!=dmpsi[2]) &&
        (dmjac[0]*dmjac[1]*dmjac[2]>1))
        mexErrMsgTxt("Deformation and Jacobian fields should have the same size.");
    if (dmjac[3]!=3 || dmjac[4]!=3)
        mexErrMsgTxt("Incompatible dimensions.");
    
    /* Dimensions of output volumes */
    if (nrhs>=4)
    {
        if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) ||
              mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[3])!= 3)
            mexErrMsgTxt("Output dimensions must have three elements");
        dm0[0] = (mwSize)floor((double)mxGetPr(prhs[3])[0]);
        dm0[1] = (mwSize)floor((double)mxGetPr(prhs[3])[1]);
        dm0[2] = (mwSize)floor((double)mxGetPr(prhs[3])[2]);
    }
    else
    {
        dm0[0] = dm1[0];
        dm0[1] = dm1[1];
        dm0[2] = dm1[2];
    }
    dm0[3]  = dm1[3];
    plhs[0] = mxCreateNumericArray(4,dm0, mxSINGLE_CLASS, mxREAL);

 
    /* Create count volume if required */
    if (nlhs>=2)
    {
        plhs[1] = mxCreateNumericArray(3,dm0, mxSINGLE_CLASS, mxREAL);
        S0      = (float *)mxGetPr(plhs[1]);
    }
    else
        S0      = (float *)0;
 
    if ((flag&2)==2)
        jpush (dm0, dm1[0]*dm1[1]*dm1[2], dm0[3], (float *)mxGetPr(prhs[1]), (float *)mxGetPr(prhs[0]), (float *)mxGetPr(plhs[0]), S0, (float *)mxGetPr(prhs[2]), dmjac[0]*dmjac[1]*dmjac[2]);
    else
        jpushc(dm0, dm1[0]*dm1[1]*dm1[2], dm0[3], (float *)mxGetPr(prhs[1]), (float *)mxGetPr(prhs[0]), (float *)mxGetPr(plhs[0]), S0, (float *)mxGetPr(prhs[2]), dmjac[0]*dmjac[1]*dmjac[2]);

}

//-------------------------------------------------------------------------
// Main MEX
//-------------------------------------------------------------------------

#include<string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    set_bound(get_bound());

    if ((nrhs>=1) && mxIsChar(prhs[0]))
    {
        int buflen;
        char *fnc_str;
        buflen = mxGetNumberOfElements(prhs[0]);
        fnc_str = (char *)mxCalloc(buflen+1,sizeof(mxChar));
        mxGetString(prhs[0],fnc_str,buflen+1);

        if (!strcmp(fnc_str,"pull"))
        {
            mxFree(fnc_str);
            jpull_mexFunction(2, nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"pullc"))
        {
            mxFree(fnc_str);
            jpull_mexFunction(0, nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"push"))
        {
            mxFree(fnc_str);
            jpush_mexFunction(2, nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"pushc"))
        {
            mxFree(fnc_str);
            jpush_mexFunction(0, nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else
        {
            mxFree(fnc_str);
            mexErrMsgTxt("Option not recognised.");
        }
    }
    else
    {
        mexErrMsgTxt("Option not recognised.");
    }
}
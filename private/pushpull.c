#include <math.h>
#include "mex.h"
#include "shoot_boundary.h"

//-------------------------------------------------------------------------
// Helpers
//-------------------------------------------------------------------------

static float min(float a, float b)
{
    if (a <= b)
        return a;
    else
        return b;
}

static float max(float a, float b)
{
    if (a >= b)
        return a;
    else
        return b;
}

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
 * //@param win   [uint 3]            Selection window (0=gauss|1=rect|2=sinc)
 */
#define TINY 5e-2f
static void jpushpull(mwSize dm0[], mwSize m1, mwSize n, 
                      float Psi[], float F0[], float S0[], float F1[], 
                      unsigned int code, float Jac[], mwSize mJ)
{
    mwSize i,                   /* Index pulled voxels */
           k,                   /* Index features */
           m0;                  /* Number of reference voxels */
    float  *px, *py, *pz;       /* Deformation components x/y/z */
    float  *jxx, *jyy, *jzz,    /* Jacobian components (diag) */ 
           *jxy, *jxz, *jyz,    /*                     (off-diag) */
           *jyx, *jzx, *jzy;    /*                     (off-diag) */
    float  NaN = mxGetNaN();

    px  = Psi;                  /* Deformation in x */
    py  = Psi+m1;               /*                y */
    pz  = Psi+m1*2;             /*                z */
    jxx = Jac;                  /* Jacobian in xx */
    jyx = Jac+mJ;               /*             yx */
    jzx = Jac+mJ*2;             /*             zx */
    jxy = Jac+mJ*3;             /*             xy */
    jyy = Jac+mJ*4;             /*             yy */
    jzy = Jac+mJ*5;             /*             zy */
    jxz = Jac+mJ*6;             /*             xz */
    jyz = Jac+mJ*7;             /*             yz */
    jzz = Jac+mJ*8;             /*             zz */
    m0  = dm0[0]*dm0[1]*dm0[2]; /* Nb of reference voxels */

    mwSize oy, oz; /* Offsets between refeence points/lines/slices */
    oy = dm0[0];
    oz = dm0[0]*dm0[1];
    
    
    float  Jxx, Jyy, Jzz, Jxy, Jxz, Jyz, Jyx, Jzx, Jzy; /* Jacobian */
    float  Cxx, Cyy, Czz, Cxy, Cxz, Cyz, Cyx, Czx, Czy; /* Covariance J*J'+I */
    float  Txx, Tyy, Tzz, Txy, Txz, Tyz, Tyx, Tzx, Tzy; /* Inverse covariance */
    float  idt;                                         /* Inverse of determinant */
    float  limx, limy, limz;                            /* Bounding box in ref space */
    float  sig    = 1./sqrt(8.*log(2.));                /* Gaussian sd so that FWHM = 1 px */
    float  pinorm = 1./sqrt(8*M_PI*M_PI*M_PI);          /* Part of kernel normalisation */
    float  norm;                                        /* Kernel normalisation */
    
    int voxJ = (mJ != 1);
    if (!voxJ)
    /* Invert Jacobian if same for all voxels */
    {
        /* Jacobian at point i */
        Jxx = *(jxx);
        Jyy = *(jyy);
        Jzz = *(jzz);
        Jxy = *(jxy);
        Jxz = *(jxz);
        Jyz = *(jyz);
        Jyx = *(jyx);
        Jzx = *(jzx);
        Jzy = *(jzy);
        
        /* Covariance = J*J'+I */
        Cxx = Jxx*Jxx + Jxy*Jxy + Jxz*Jxz + 1;
        Cxy = Jxx*Jyx + Jxy*Jyy + Jxz*Jyz;
        Cxz = Jxx*Jzx + Jxy*Jzy + Jxz*Jzz;
        Cyx = Jyx*Jxx + Jyy*Jxy + Jyz*Jxz;
        Cyy = Jyx*Jyx + Jyy*Jyy + Jyz*Jyz + 1;
        Cyz = Jyx*Jzx + Jyy*Jzy + Jyz*Jzz;
        Czx = Jzx*Jxx + Jzy*Jxy + Jzz*Jxz;
        Czy = Jzx*Jyx + Jzy*Jyy + Jzz*Jyz;
        Czz = Jzx*Jzx + Jzy*Jzy + Jzz*Jzz + 1;
        
        /* Invert covariance */
        idt = 1.0/(Cxx*Cyy*Czz + 2*Cxy*Cxz*Cyz
                   -Cxx*Cyz*Cyz-Cyy*Cxz*Cxz-Czz*Cxy*Cxy);
        norm = sqrt(idt)*pinorm/(sig*sig);
        Txx = idt*(Cyy*Czz-Cyz*Cyz)/(sig*sig);
        Tyx = idt*(Cxz*Cyz-Cxy*Czz)/(sig*sig);
        Tzx = idt*(Cxy*Cyz-Cxz*Cyy)/(sig*sig);
        Txy = idt*(Cxz*Cyz-Cxy*Czz)/(sig*sig);
        Tyy = idt*(Cxx*Czz-Cxz*Cxz)/(sig*sig);
        Tzy = idt*(Cxy*Cxz-Cxx*Cyz)/(sig*sig);
        Txz = idt*(Cxy*Cyz-Cxz*Cyy)/(sig*sig);
        Tyz = idt*(Cxy*Cxz-Cxx*Cyz)/(sig*sig);
        Tzz = idt*(Cxx*Cyy-Cxy*Cxy)/(sig*sig);
        
        /* Bounding box of contributing points in ref space */
        limx = limy = limz = 0;
        limx = max(limx, fabs(Jxx));
        limx = max(limx, fabs(Jxy));
        limx = max(limx, fabs(Jxz));
        limy = max(limy, fabs(Jyx));
        limy = max(limy, fabs(Jyy));
        limy = max(limy, fabs(Jyz));
        limz = max(limz, fabs(Jzx));
        limz = max(limz, fabs(Jzy));
        limz = max(limz, fabs(Jzz));
        limx += 1;
        limy += 1;
        limz += 1;
        limx *= 3*sig;
        limy *= 3*sig;
        limz *= 3*sig;
        mexPrintf("Limits: %f %f %f\n", limx, limy, limz);
    }
    
    float *pf1 = F1;
    for (i=0; i<m1; ++i, ++pf1)    /* loop over pulled voxels */
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
                /* Jacobian at point i */
                Jxx = *(jxx++);
                Jyy = *(jyy++);
                Jzz = *(jzz++);
                Jxy = *(jxy++);
                Jxz = *(jxz++);
                Jyz = *(jyz++);
                Jyx = *(jyx++);
                Jzx = *(jzx++);
                Jzy = *(jzy++);
                
                /* Covariance = J*J'+I */
                Cxx = Jxx*Jxx + Jxy*Jxy + Jxz*Jxz + 1;
                Cxy = Jxx*Jyx + Jxy*Jyy + Jxz*Jyz;
                Cxz = Jxx*Jzx + Jxy*Jzy + Jxz*Jzz;
                Cyx = Jyx*Jxx + Jyy*Jxy + Jyz*Jxz;
                Cyy = Jyx*Jyx + Jyy*Jyy + Jyz*Jyz + 1;
                Cyz = Jyx*Jzx + Jyy*Jzy + Jyz*Jzz;
                Czx = Jzx*Jxx + Jzy*Jxy + Jzz*Jxz;
                Czy = Jzx*Jyx + Jzy*Jyy + Jzz*Jyz;
                Czz = Jzx*Jzx + Jzy*Jzy + Jzz*Jzz + 1;

                /* Invert covariance */
                idt = 1.0/(Cxx*Cyy*Czz + 2*Cxy*Cxz*Cyz
                           -Cxx*Cyz*Cyz-Cyy*Cxz*Cxz-Czz*Cxy*Cxy);
                norm = sqrt(idt)*pinorm/(sig*sig);
                Txx = idt*(Cyy*Czz-Cyz*Cyz)/(sig*sig);
                Tyx = idt*(Cxz*Cyz-Cxy*Czz)/(sig*sig);
                Tzx = idt*(Cxy*Cyz-Cxz*Cyy)/(sig*sig);
                Txy = idt*(Cxz*Cyz-Cxy*Czz)/(sig*sig);
                Tyy = idt*(Cxx*Czz-Cxz*Cxz)/(sig*sig);
                Tzy = idt*(Cxy*Cxz-Cxx*Cyz)/(sig*sig);
                Txz = idt*(Cxy*Cyz-Cxz*Cyy)/(sig*sig);
                Tyz = idt*(Cxy*Cxz-Cxx*Cyz)/(sig*sig);
                Tzz = idt*(Cxx*Cyy-Cxy*Cxy)/(sig*sig);

                /* Bounding box of contributing points in ref space */
                limx = limy = limz = 0;
                limx = max(limx, fabs(Jxx));
                limx = max(limx, fabs(Jxy));
                limx = max(limx, fabs(Jxz));
                limy = max(limy, fabs(Jyx));
                limy = max(limy, fabs(Jyy));
                limy = max(limy, fabs(Jyz));
                limz = max(limz, fabs(Jzx));
                limz = max(limz, fabs(Jzy));
                limz = max(limz, fabs(Jzz));
                limx += 1;
                limy += 1;
                limz += 1;
                limx *= 3*sig;
                limy *= 3*sig;
                limz *= 3*sig;
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
            for(k=0; k<lx; ++k) oox[k] = bound(ixmin+k, dm0[0]);
            for(k=0; k<ly; ++k) ooy[k] = bound(iymin+k, dm0[1])*oy;
            for(k=0; k<lz; ++k) ooz[k] = bound(izmin+k, dm0[2])*oz;
            float  ddx[128], ddy[128], ddz[128];
            for(k=0; k<lx; ++k) ddx[k] = x-(float)(k+ixmin);
            for(k=0; k<ly; ++k) ddy[k] = y-(float)(k+iymin);
            for(k=0; k<lz; ++k) ddz[k] = z-(float)(k+izmin);
            
            mwSignedIndex  ix, iy, iz;          /* Reference voxel indices */
            float   dx, dy, dz;                 /* Distance to pulled position */
            float   Tdx, Tdy, Tdz;              /* Transformed distance */
            float   Td2;                        /* Squared distance */
            double  acc[n];                     /* Accumulator per feature */
            float   *pf0x, *pf0y, *pf0z;        /* Pointer into ref volume */
            float   *ps0x, *ps0y, *ps0z;        /* Pointer into count volume */
            float   *pf1k;                      /* Pointer into subject volume */
            double  *pacc;
            float   w;
            
            if ((code&1)==1)
            /* Pull */
            {
                /* Initialise accumulator */
                pacc = acc;
                for (k=0; k<n; ++k)
                    *(pacc++) = 0.;
                            
                /* Loop over contributing voxels in ref space */
                for (iz=0; iz<lz;)
                {
                    dz   = ddz[iz];
                    pf0z = F0 + ooz[iz++];
                    for (iy=0; iy<ly;)
                    {
                        dy   = ddy[iy];
                        pf0y = pf0z + ooy[iy++];
                        for (ix=0; ix<lx;)
                        {
                            dx   = ddx[ix];
                            pf0x = pf0y + oox[ix++];
                            
                            Tdx  = Txx * dx + Txy * dy + Txz * dz;
                            Tdy  = Tyx * dx + Tyy * dy + Tyz * dz;
                            Tdz  = Tzx * dx + Tzy * dy + Tzz * dz;
                            Td2  = dx * Tdx + dy * Tdy + dz * Tdz;
                            if (Td2 < 9.)
                            {
                                w    = exp(-Td2)*norm;
                            
                                pacc = acc;
                                for (k=0; k<n; ++k, pf0x += m0) {
                                    if (pf0x-F0 <0 || pf0x-F0 >= m0*n)
                                        mexErrMsgTxt("Out of bound! (pull acc)");
                                    *(pacc++) += w * (*pf0x);
                                }
                            }
                        }
                    }
                }

                /* Fill output volume */
                pf1k = pf1;
                pacc = acc;
                for (k=0; k<n; ++k, pf1k += m1) {
                    if (pf1k-F1 <0 || pf1k-F1 >= m1*n)
                        mexErrMsgTxt("Out of bound! (pull fill).");
                    (*pf1k) = *(pacc++);
                }
            }
            else
            /* Push */
            {
                /* Loop over contributing voxels in ref space */
                for (iz=0; iz<lz;)
                {
                    dz   = z-(float)(iz+izmin);
                    pf0z = F0 + ooz[iz];
                    ps0z = S0 + ooz[iz++];
                    for (iy=0; iy<ly;)
                    {
                        dy   = y-(float)(iy+iymin);
                        pf0y = pf0z + ooy[iy];
                        ps0y = ps0z + ooy[iy++];
                        for (ix=0; ix<lx;)
                        {
                            dx   = x-(float)(ix+ixmin);
                            pf0x = pf0y + oox[ix];
                            ps0x = ps0y + oox[ix++];
                            
                            Tdx  = Txx * dx + Txy * dy + Txz * dz;
                            Tdy  = Tyx * dx + Tyy * dy + Tyz * dz;
                            Tdz  = Tzx * dx + Tzy * dy + Tzz * dz;
                            Td2  = dx * Tdx + dy * Tdy + dz * Tdz;
                            if (Td2 < 9.)
                            {
                                w    = exp(-Td2)*norm;

                                pf1k = pf1;
                                for (k=0; k<n; ++k, pf0x += m0, pf1k += m1) {
                                    if (pf0x-F0 <0 || pf0x-F0 >= m0*n)
                                        mexErrMsgTxt("Out of bound! (push F0).");
                                    if (pf1k-F1 <0 || pf1k-F1 >= m1*n)
                                        mexErrMsgTxt("Out of bound! (push F1).");
                                    (*pf0x) += w * (*pf1k);
                                }
                                /* Count image */
                                if (S0!=0) {
                                    if (ps0x-S0 <0 || ps0x-S0 >= m0)
                                        mexErrMsgTxt("Out of bound! (push S0).");
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
            for(k=0; k<n; ++k, pf1k += m1)
            {
                if (pf1k-F1 <0 || pf1k-F1 >= m1*n)
                    mexErrMsgTxt("Out of bound! (NaN).");
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
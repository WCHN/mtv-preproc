function nfname = downsample_inplane(fname,vxd)
% Down-sample a NIfTI image in the high-resolution plane
% FORMAT nfname = downsample_inplane(fname)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging  
% Get image data

if nargin < 2, vxd = [1 1 1]; end

Nii = nifti(fname);
M0  = Nii.mat;             
X   = Nii.dat(:,:,:);  
dm0 = size(X);                                                   
vx0 = sqrt(sum(M0(1:3,1:3).^2));

% Get down-sampling factor
d = vx0(1:2);
if d(1)>=1, d(1) = vxd(1); end
if d(2)>=1, d(2) = vxd(2); end
d(3) = 1;

if round(d(1),3)<1 || round(d(2),3)<1        
    % NN downsampling
    d = d./vxd;
    D = diag([d, 1]);          
           
    dm1 = floor(D(1:3,1:3)*dm0')';
    M1  = M0/D;       
    
    T = M0\M1;
    y = make_deformation(T,dm1);
                   
    X = spm_bsplins(X,y(:,:,:,1),y(:,:,:,2),y(:,:,:,3),[0 0 0 0 0 0]);
    clear y

    X(~isfinite(X)) = 0;
else
    M1 = M0;
end

fname         = Nii.dat.fname;
[pth,nam,ext] = fileparts(fname);
nfname        = fullfile(pth,['ds' nam ext]);
        
create_nii(nfname,X,M1,Nii.dat.dtype,Nii.descrip,Nii.dat.offset,Nii.dat.scl_slope,Nii.dat.scl_inter);
%==========================================================================

%==========================================================================
function y = make_deformation(M,dm)
[x0,y0,z0] = ndgrid(1:dm(1),...
                    1:dm(2),...
                    1:dm(3));
y          = cat(4,M(1,1)*x0 + M(1,2)*y0 + M(1,3)*z0 + M(1,4), ...
                   M(2,1)*x0 + M(2,2)*y0 + M(2,3)*z0 + M(2,4), ...
                   M(3,1)*x0 + M(3,2)*y0 + M(3,3)*z0 + M(3,4));
%==========================================================================  
function [img,mat,dm] = resample_img(Nii,samp,deg,bc)
% Resample an image using deg interpolation, with bc boundary conditions.
% If samp < 1, does down-sampling; if samp > 1, does up-sampling.
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, deg = 0; end
if nargin < 4, bc  = 0; end

if numel(samp) == 1, samp = samp*ones([1 3]); end
if numel(deg)  == 1, deg  = deg*ones([1 3]);  end
if numel(bc)   == 1, bc   = bc*ones([1 3]);   end

if numel(Nii.dat.dim) ~= 3
    error('numel(Nii.dat.dim) ~= 3')
end

% Input image properties
img  = Nii.dat(:,:,:);    
mat0 = Nii.mat;
dm0  = size(img);

% Output image properties
D    = diag([samp 1]);
mat  = mat0/D;
dm   = floor(D(1:3,1:3)*dm0')';

% Make interpolation grid
[x0,y0,z0] = ndgrid(1:dm(1),1:dm(2),1:dm(3));

T = mat0\mat;    

x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

% Resample
img                         = spm_bsplins(img,x1,y1,z1,[deg bc]);    
img(~isfinite(img) | img<0) = 0;
%==========================================================================
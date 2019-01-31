function [y,msk] = get_y_superres(Nii,dat,dm,mat,deg)
% Get super-resolution starting estimates using bspline interpolation
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 5, deg = 4; end

N          = numel(Nii);
bs         = [deg deg deg  0 0 0];
[x0,y0,z0] = ndgrid(1:dm(1),1:dm(2),1:dm(3));

y = single(0);
for n=1:N

    T = dat.A(n).mat\mat;    

    x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
    y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
    z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

    coeff = spm_bsplinc(get_nii(Nii(n)),bs);    
    img   = spm_bsplins(coeff,x1,y1,z1,bs);
    
    img(~isfinite(img) | img<0) = 0;

    y = y + single(img); 
end
y   = y/N;
msk = isfinite(y) & y ~= 0;
%==========================================================================
function [psnr1,ssim1] = compute_image_metrics(Nii,Nii_ref)
% Compute image quality metrics
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimagin

C = numel(Nii_ref);

if iscell(Nii(1))
    img = get_nii(Nii{1});
else
    img = get_nii(Nii(1));
end
ref = get_nii(Nii_ref(1));      

for c=2:C
    if iscell(Nii(c))
        img1 = get_nii(Nii{c});
    else
        img1 = get_nii(Nii(c));
    end

    img = cat(3,img,img1);
    ref = cat(3,ref,get_nii(Nii_ref(c)));                    
end

psnr1 = get_psnr(img(:),ref(:));
ssim1 = ssim(img,ref);
%==========================================================================
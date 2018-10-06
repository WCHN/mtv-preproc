function put_nii(nii,c,img)
% Write to NIfTI
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

dm = nii(c).dat.dim;
if numel(dm) == 4
    nii(c).dat(:,:,:,:) = img;
elseif numel(dm) == 3
    nii(c).dat(:,:,:) = img;
end
%==========================================================================
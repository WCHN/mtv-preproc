function put_nii(nii,img)
% Write to NIfTI
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

dm = nii.dat.dim;
if numel(dm) == 4
    nii.dat(:,:,:,:) = img;
elseif numel(dm) == 3
    nii.dat(:,:,:) = img;
end
%==========================================================================
function nii = put_nii(nii,img)
% Write data
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

dm = size(img);
if numel(dm) == 4
    nii.dat(:,:,:,:) = img;
elseif numel(dm) == 3
    nii.dat(:,:,:) = img;
end
%==========================================================================
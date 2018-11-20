function img = get_nii(nii,z)
% Read data
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 2, z = 0; end

if z == 0
    img = single(nii.dat(:,:,:,:,:));
else
    img = single(nii.dat(:,:,z,:,:));
end
%==========================================================================
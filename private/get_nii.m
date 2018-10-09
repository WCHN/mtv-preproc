function img = get_nii(nii,z)
% Read from NIfTI
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 2, z = 0; end

dm = nii.dat.dim;
if z == 0
    if numel(dm) == 4
        img = single(nii.dat(:,:,:,:));
    elseif numel(dm) == 3
        img = single(nii.dat(:,:,:));
    end
else
    if numel(dm) == 4
        img = single(nii.dat(:,:,z,:));
    elseif numel(dm) == 3
        img = single(nii.dat(:,:,z));
    end
end
%==========================================================================
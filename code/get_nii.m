function img = get_nii(nii,c,z)
% Read from NIfTI
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, z = 0; end

dm = nii(c).dat.dim;
if z == 0
    if numel(dm) == 4
        img = single(nii(c).dat(:,:,:,:));
    elseif numel(dm) == 3
        img = single(nii(c).dat(:,:,:));
    end
else
    if numel(dm) == 4
        img = single(nii(c).dat(:,:,z,:));
    elseif numel(dm) == 3
        img = single(nii(c).dat(:,:,z));
    end
end
%==========================================================================
function img = get_nii(nii,z,dim)
% Read data
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 2, z   = 0; end
if nargin < 3, dim = 3; end

if z == 0
    img = single(nii.dat());
else
    img = single(select_slices(nii.dat,dim,z));
end
%==========================================================================
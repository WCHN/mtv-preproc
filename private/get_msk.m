function msk = get_msk(f,mu)
% Get mask defining what voxels to regard as missing data
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

msk = isfinite(f) & isfinite(mu) & f ~= 0;
msk = msk(:);
%==========================================================================
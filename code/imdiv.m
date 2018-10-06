function out = imdiv(img,lam,dm,vx)
% Compute image divergence
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% out = single(full(reshape(lam*D'*double(img),dm)));
out = lam*spm_imbasics('dive',img(:,:,:,1),img(:,:,:,2),img(:,:,:,3),vx);
%==========================================================================
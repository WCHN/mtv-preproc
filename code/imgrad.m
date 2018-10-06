function out = imgrad(img,lam,dm,vx)
% Compute image gradient
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% out = single(full(reshape(lam*D*double(img(:)),[dm 3])));
out = zeros([dm 3],'single');
[out(:,:,:,1),out(:,:,:,2),out(:,:,:,3)] = spm_imbasics('grad',img,vx);
out = lam*out;
%==========================================================================
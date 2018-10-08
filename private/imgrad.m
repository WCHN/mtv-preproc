function out = imgrad(img,lam,dm,vx)
% Compute image gradient
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% out = single(full(reshape(lam*D*double(img(:)),[dm 3])));
out = zeros([dm 3],'single');
[out(:,:,:,1),out(:,:,:,2),out(:,:,:,3)] = grad(img,vx);
out = lam*out;
%==========================================================================

%==========================================================================
function [Dx,Dy,Dz] = grad(X,vx) 
% Calculate 2D or 3D gradient of an image (with voxel size)
% FORMAT [Dx,Dy,Dz] = grad(X,vx)
% X          - Image
% vx         - voxel size
% [Dx,Dy,Dz] - Gradients in x-,y- and z-direction
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   
if nargin < 2, vx = ones([1 3],'like',X); end

if size(X,3)==1
    Dx = [diff(X,1,2),zeros(size(X,1),1,'like',X)]./vx(1);
    Dy = [diff(X,1,1);zeros(1,size(X,2),'like',X)]./vx(2);
    Dz = 0;
else
    Dx = cat(2,diff(X,1,2),zeros(size(X,1),1,size(X,3),'like',X))./vx(1);
    Dy = cat(1,diff(X,1,1),zeros(1,size(X,2),size(X,3),'like',X))./vx(2);
    Dz = cat(3,diff(X,1,3),zeros(size(X,1),size(X,2),1,'like',X))./vx(3);  
end
%==========================================================================
function out = imdiv(img,lam,dm,vx)
% Compute image divergence
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% out = single(full(reshape(lam*D'*double(img),dm)));
out = lam*dive(img(:,:,:,1),img(:,:,:,2),img(:,:,:,3),vx);
%==========================================================================

%==========================================================================
function div = dive(Dx,Dy,Dz,vx)  
% Computes the divergence of an image (with voxel size)
% FORMAT div = dive(Dx,Dy,Dz,vx) 
% [Dx,Dy,Dz] - Gradients in x-,y- and z-direction
% vx         - Voxel size
% div        - Divergence
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin < 4, vx = ones([1 3],'like',Dx); end

if size(Dx,3) == 1
    Du  = [-Dx(:,1), -diff(Dx(:,1:end-1),1,2), Dx(:,end-1)];
    Dv  = [-Dy(1,:); -diff(Dy(1:end-1,:),1,1); Dy(end-1,:)];
    div = Du./vx(1) + Dv./vx(2);
else
    Du  = cat(2, -Dx(:,1,:), -diff(Dx(:,1:end-1,:),1,2), Dx(:,end-1,:)); 
    Dv  = cat(1, -Dy(1,:,:), -diff(Dy(1:end-1,:,:),1,1), Dy(end-1,:,:));
    Dw  = cat(3, -Dz(:,:,1), -diff(Dz(:,:,1:end-1),1,3), Dz(:,:,end-1));
    div = Du./vx(1) + Dv./vx(2) + Dw./vx(3);
end
%==========================================================================
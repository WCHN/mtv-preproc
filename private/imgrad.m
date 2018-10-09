function G = imgrad(im,vx) 
% Compute image gradient
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   

if nargin < 2, vx = ones([1 3],'like',im); end

Gx = cat(2, diff(im,1,2), zeros(size(im,1),1,size(im,3), 'like', im))./vx(1);
Gy = cat(1, diff(im,1,1), zeros(1,size(im,2),size(im,3), 'like', im))./vx(2);
Gz = cat(3, diff(im,1,3), zeros(size(im,1),size(im,2),1, 'like', im))./vx(3);  

G = cat(4,Gx,Gy,Gz);
%==========================================================================
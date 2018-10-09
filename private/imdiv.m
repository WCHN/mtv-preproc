function D = imdiv(imx,imy,imz,vx)  
% Compute image divergence
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 4, vx = ones([1 3],'like',imx); end

Du = cat(2, -imx(:,1,:), -diff(imx(:,1:end-1,:),1,2), imx(:,end-1,:)); 
Dv = cat(1, -imy(1,:,:), -diff(imy(1:end-1,:,:),1,1), imy(end-1,:,:));
Dw = cat(3, -imz(:,:,1), -diff(imz(:,:,1:end-1),1,3), imz(:,:,end-1));

D = Du./vx(1) + Dv./vx(2) + Dw./vx(3);
%==========================================================================
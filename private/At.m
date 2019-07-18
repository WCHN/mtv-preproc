function Y = At(X,dat,tau,n)  
% Adjoint of forward model (y=A'x)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

USE_PUSHPULL = true;

if nargin < 3, tau = ones(1,dat.N); end
    
if nargin < 4
    Y = single(0);   
    for n=1:dat.N      
%         R = dat.A(n).R;
%         T = dat.mat\R*dat.A(n).mat;
%         y = apply_affine(T,dat.A(n).dm);
% 
%         if USE_PUSHPULL
%             % Using pushpull.c
%             if strcmp(dat.method,'superres')
%                 tmp = pushpull('pushc',X{n},y,single(dat.A(n).J),double(dat.A(n).win),double(dat.dm));     
%             elseif strcmp(dat.method,'denoise')
%                 tmp = spm_diffeo('pushc',X{n},y,dat.dm);
%             end
%             clear y
%             tmp(~isfinite(tmp)) = 0;
%         else
%             % Using spm_diffeo
%             tmp                 = spm_diffeo('push',X{n},y,dat.dm);     
%             clear y
%             tmp(~isfinite(tmp)) = 0;    
%             tmp                 = real(ifftn(fftn(tmp).*dat.A(n).S)); 
%             tmp(~isfinite(tmp)) = 0;
%         end
        tmp = apply_proj(X{n},dat,n,'At');
        Y   = Y + tau(n).*tmp;           
    end
else
%     R = dat.A(n).R;
%     T = dat.mat\R*dat.A(n).mat;
%     y = apply_affine(T,dat.A(n).dm);
% 
%     if USE_PUSHPULL
%         % Using pushpull.c
%         if strcmp(dat.method,'superres')
%             tmp = pushpull('pushc',X,y,single(dat.A(n).J),double(dat.A(n).win),double(dat.dm));     
%         elseif strcmp(dat.method,'denoise')
%             tmp = spm_diffeo('pushc',X,y,dat.dm);
%         end
%         clear y
%         tmp(~isfinite(tmp)) = 0;
%     else
%         % Using spm_diffeo
%         tmp                 = spm_diffeo('push',X,y,dat.dm);     
%         clear y
%         tmp(~isfinite(tmp)) = 0;    
%         tmp                 = real(ifftn(fftn(tmp).*dat.A(n).S)); 
%         tmp(~isfinite(tmp)) = 0;            
%     end
    
    tmp = apply_proj(X,dat,n,'At');
    Y   = tau(n).*tmp;      
end
%==========================================================================
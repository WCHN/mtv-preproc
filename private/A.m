function X = A(Y,dat,n)
% Forward model (x=Ay)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

USE_PUSHPULL = true;

if ~USE_PUSHPULL
    Y(~isfinite(Y)) = 0;
    Y               = fftn(Y);
end

if nargin < 3
    X = cell(1,dat.N);
    for n=1:dat.N
%         X{n} = apply_proj(Y,dat,n,'A');
%         R = dat.A(n).R;
%         T = dat.mat\R*dat.A(n).mat;
%         y = apply_affine(T,dat.A(n).dm);
% 
%         if USE_PUSHPULL
%             % Using pushpull.c
%             if strcmp(dat.method,'superres')
%                 X{n} = pushpull('pullc',Y,y,single(dat.A(n).J),double(dat.A(n).win));    
%             elseif strcmp(dat.method,'denoise')
%                 X{n} = spm_diffeo('pullc',Y,y);
%             end
%             clear y    
%             X{n}(~isfinite(X{n})) = 0; 
%         else
%             % Using spm_diffeo
%             tmp                   = real(ifftn(Y.*dat.A(n).S));    
%             X{n}                  = spm_diffeo('pull',tmp,y);    
%             X{n}(~isfinite(X{n})) = 0; 
%             clear y tmp
%         end

        X{n} = apply_proj(Y,dat,n,'A');
    end
else
%     R = dat.A(n).R;
%     T = dat.mat\R*dat.A(n).mat;
%     y = apply_affine(T,dat.A(n).dm);
% 
%     if USE_PUSHPULL
%         % Using pushpull.c
%         if strcmp(dat.method,'superres')
%             X = pushpull('pullc',Y,y,single(dat.A(n).J),double(dat.A(n).win));    
%         elseif strcmp(dat.method,'denoise')
%             X = spm_diffeo('pullc',Y,y);
%         end
%         clear y    
%     else
%         % Using spm_diffeo
%         tmp             = real(ifftn(Y.*dat.A(n).S));    
%         X               = spm_diffeo('pull',tmp,y);    
%         X(~isfinite(X)) = 0; 
%         clear y    
%     end

    X  = apply_proj(Y,dat,n,'A');
%         Y1 = apply_proj(X,dat,n,'At');
end
%==========================================================================  
function Y = At(X,dat,tau)  
% Adjoint of forward model (y=A'x)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, tau = ones(1,dat.N); end
    
Y = single(0);   
for n=1:dat.N      
    R = dat.A(n).R;
    T = dat.mat\R*dat.A(n).mat;
    y = apply_affine(T,dat.A(n).dm);
    
    if strcmp(dat.method,'superres')
        tmp = pushpull('push',X{n},y,single(dat.A(n).J),double(dat.A(n).win),double(dat.dm));     
    elseif strcmp(dat.method,'denoise')
        tmp = spm_diffeo('push',X{n},y,dat.dm);
    end
    clear y
     
    Y = Y + tau(n).*tmp;           
end
%==========================================================================
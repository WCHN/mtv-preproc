function Y1 = AtA(Y,prior,tau,lam,dat)  
% A'A
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

F               = Y;        
F(~isfinite(F)) = 0;
F               = fftn(F);    

Y1 = single(0);
for n=1:dat.N
    T   = dat.mat\dat.A(n).mat;
    y   = apply_affine(T,dat.A(n).dm);
        
    tmp = real(ifftn(F.*dat.A(n).S));    
    
    tmp                 = spm_diffeo('pull',tmp,y);
    tmp(~isfinite(tmp)) = 0; 
    
    vx1 = sqrt(sum(dat.mat(1:3,1:3).^2));
    vx0 = sqrt(sum(dat.A(n).mat(1:3,1:3).^2));
    scl = prod(vx1./vx0);
    tmp = scl*tmp;
    
    tmp                 = spm_diffeo('push',tmp,y,dat.dm);          
    clear y
    tmp(~isfinite(tmp)) = 0;    
    
    tmp = real(ifftn(fftn(tmp).*dat.A(n).S));
    
    Y1  = Y1 + tau(n).*tmp;
end

Y1 = Y1 + lam*prior(Y);
%==========================================================================
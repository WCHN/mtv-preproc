function Y1 = AtA(Y,prior,tau,lam,dat)  
% A'A
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
  
Y(~isfinite(Y)) = 0;   

Y1 = single(0);
for n=1:dat.N
    T   = dat.mat\dat.A(n).mat;
    y   = apply_affine(T,dat.A(n).dm);
        
    tmp = Y;
    
    tmp                 = pushpull('pull',single(tmp),single(y),single(dat.A(n).J));
    tmp(~isfinite(tmp)) = 0; 
    
    tmp                 = pushpull('push',single(tmp),single(y),single(dat.A(n).J),double(dat.dm));          
    clear y
    tmp(~isfinite(tmp)) = 0;    
    
    Y1  = Y1 + tau(n).*tmp;
end

Y1 = Y1 + lam*prior(Y);
%==========================================================================
function Y = At(X,dat,tau)  
% Adjoint of forward model (y=A'x)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, tau = ones(1,dat.N); end
    
Y = single(0);   
for n=1:dat.N      
    T = dat.mat\dat.A(n).mat;
    y = apply_affine(T,dat.A(n).dm);
    
    tmp = spm_diffeo('push',X{n},y,dat.dm);     
    clear y
    tmp(~isfinite(tmp)) = 0;
       
    tmp = real(ifftn(fftn(tmp).*dat.A(n).S));  
     
    Y = Y + tau(n).*tmp;           
end
%==========================================================================
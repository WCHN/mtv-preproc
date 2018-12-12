function Y = At(X,dat,tau)  
% Adjoint of forward model (y=A'x)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, tau = ones(1,dat.N); end
    
Y = single(0);   
for n=1:dat.N      
    T = dat.mat\dat.A(n).mat;
    y = apply_affine(T,dat.A(n).dm);
    
    tmp = pushpull('push',single(X{n}),single(y),single(dat.A(n).J),double(dat.dm));     
    clear y
    tmp(~isfinite(tmp)) = 0;
     
    Y = Y + tau(n).*tmp;           
end
%==========================================================================
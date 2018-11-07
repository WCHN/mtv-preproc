function X = A(Y,dat)
% Forward model (x=Ay)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

Y(~isfinite(Y)) = 0;
Y               = fftn(Y);    

X = cell(1,dat.N);
for n=1:dat.N
    T    = dat.mat\dat.A(n).mat;
    y    = apply_affine(T,dat.A(n).dm);
     
    tmp  = real(ifftn(Y.*dat.A(n).S));  
    
    X{n}                  = spm_diffeo('pull',tmp,y);    
    X{n}(~isfinite(X{n})) = 0; 
    clear y tmp
    
    X{n} = X{n};
end
%==========================================================================  
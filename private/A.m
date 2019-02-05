function X = A(Y,dat)
% Forward model (x=Ay)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

Y(~isfinite(Y)) = 0; 

X = cell(1,dat.N);
for n=1:dat.N
    T    = dat.mat\dat.A(n).mat;
    y    = apply_affine(T,dat.A(n).dm);
    
    X{n}                  = pushpull('pull',single(Y),single(y),single(dat.A(n).J),double(dat.A(n).win));    
    X{n}(~isfinite(X{n})) = 0; 
    clear y tmp
end
%==========================================================================  
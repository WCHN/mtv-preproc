function X = A(Y,dat)
% Forward model (x=Ay)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Get rigid basis
B = get_rigid_basis;

Y(~isfinite(Y)) = 0; 

X = cell(1,dat.N);
for n=1:dat.N
    R = spm_dexpm(dat.A(n).q,B);
    T = dat.mat\R*dat.A(n).mat;
    y = apply_affine(T,dat.A(n).dm);
    
    if strcmp(dat.method,'superres')
        X{n} = pushpull('pull',single(Y),single(y),single(dat.A(n).J),double(dat.A(n).win));    
    elseif strcmp(dat.method,'denoise')
        X{n} = spm_diffeo('pull',single(Y),single(y));
    end
    clear y    
    X{n}(~isfinite(X{n})) = 0;     
end
%==========================================================================  
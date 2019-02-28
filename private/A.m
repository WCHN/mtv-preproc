function X = A(Y,dat)
% Forward model (x=Ay)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

X = cell(1,dat.N);
for n=1:dat.N
    R = dat.A(n).R;
    T = dat.mat\R*dat.A(n).mat;
    y = apply_affine(T,dat.A(n).dm);
    
    if strcmp(dat.method,'superres')
        X{n} = pushpull('pull',Y,y,single(dat.A(n).J),double(dat.A(n).win));    
    elseif strcmp(dat.method,'denoise')
        X{n} = spm_diffeo('pull',Y,y);
    end
    clear y    
end
%==========================================================================  
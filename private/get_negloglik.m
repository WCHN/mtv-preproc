function ll = get_negloglik(method,y,x,tau,dat)
% Compute negative log of likelihood part
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if strcmpi(method,'superres')    
    Ay = A(y,dat);
    ll = 0;
    for n=1:dat.N
        msk = isfinite(x{n}) & x{n} ~= 0;
        ll  = ll + 0.5*tau(n)*sum((x{n}(msk) - Ay{n}(msk)).^2);
    end
elseif strcmpi(method,'denoise')    
    msk = isfinite(x{1}) & x{1} ~= 0;
    ll  = 0.5*tau*sum(sum(sum((x{1}(msk) - y(msk)).^2)));
end
%==========================================================================
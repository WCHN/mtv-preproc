function ll = get_ll1(y,x,tau,dat)
% Compute log of likelihood part
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
  
Ay = A(y,dat);
ll = 0;
for n=1:dat.N
    msk = isfinite(x{n}) & x{n} ~= 0;
    msk = msk(:);
    ll  = ll - 0.5*tau(n)*sum((double(x{n}(msk)) - double(Ay{n}(msk))).^2);
end
%==========================================================================
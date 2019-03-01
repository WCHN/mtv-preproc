function ll = get_ll1(use_projmat,y,x,tau,dat)
% Compute log of likelihood part
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
  
if use_projmat
    % We use the projection matrices (A, At)
    
    Ay = A(y,dat);
    ll = 0;
    for n=1:dat.N
        msk = get_msk(x{n},Ay{n});
%         ll  = ll - 0.5*tau(n)*sum((double(x{n}(msk)) - double(Ay{n}(msk))).^2);
        ll  = ll - 0.5*tau(n)*sum(double(Ay{n}(msk).^2 - 2*Ay{n}(msk).*x{n}(msk)));
    end
else
    % We do not use the projection matrices (A, At)
    msk = isfinite(x{1}) & isfinite(y);
    msk = msk(:);
%     ll  = -0.5*tau*sum((double(x{1}(msk)) - double(y(msk))).^2);
    ll  = -0.5*tau*sum(double(y(msk).^2 - 2*y(msk).*x{1}(msk)));
end   
%==========================================================================
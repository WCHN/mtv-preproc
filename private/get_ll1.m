function ll1 = get_ll1(method,y,x,tau,dat,b)
% Compute likelihood term of posterior
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 6
    b = 1;
end

if strcmpi(method,'superres')    
    Ay  = A(y,dat);
    ll1 = 0;
    for n=1:numel(Ay)
        ll1 = ll1 - (tau(n)/2)*sum(sum(sum((b.*Ay{n} - x).^2)));
    end
elseif strcmpi(method,'denoise')    
    ll1 = -(tau/2)*sum(sum(sum((b.*y - x).^2)));
end
%==========================================================================
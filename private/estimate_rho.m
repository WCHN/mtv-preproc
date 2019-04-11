function rho = estimate_rho(tau,lam)
% Estimate rho (this value seems to lead to reasonably good convergence)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging    

C    = numel(tau);
atau = [];
for c=1:C
    N = numel(tau{c});
    for n=1:N
        atau = [atau tau{c}(n)];
    end
end
rho = sqrt(mean(atau))/mean(lam);        
%==========================================================================    
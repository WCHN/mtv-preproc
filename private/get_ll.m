function ll = get_ll(method,y,x,tau,dat,varargin)
% Compute likelihood term of posterior
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if numel(varargin) == 0
    % Compute negative log of likelihood part
    if strcmpi(method,'superres')    
        Ay = A(y,dat);
        ll = 0;
        for n=1:dat.N
            ll = ll + (tau(n)/2)*sum(sum(sum((Ay{n} - x{n}).^2)));
        end
    elseif strcmpi(method,'denoise')    
        ll = (tau/2)*sum(sum(sum((y - x{1}).^2)));
    end
else        
    % Parse input
    lam = varargin{1};
    u   = varargin{2};
    w   = varargin{3};
    rho = varargin{4};
    vx  = varargin{5};
    
    % Compute parts of objective function (-ln(...))
    ll = zeros(1,3);
    
    Ay  = A(y,dat);
    for n=1:numel(Ay)
        ll(1) = ll(1) + (tau(n)/(2*rho))*sum(sum(sum((Ay{n} - x{n}).^2)));
    end
    clear Ay x
    
    m     = spm_field('vel2mom',y,[vx 0 lam^2 0]);
    ll(2) = 0.5*y(:)'*m(:);
    clear m
    
    Dy    = imgrad(y,vx);
    ll(3) = lam*(w(:)/rho - u(:))'*Dy(:);
end
%==========================================================================
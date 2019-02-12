function [Nii_y,ll1,ll2,msk] = estimate_initial_y(Nii_x,Nii_y,Nii_H,dat,tau,rho,lam,vx,dm,num_workers,p)
% Compute initial estimate of recovered image(s)
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Some parameters from options struct
modality = p.Results.Modality;
method   = p.Results.Method;
nitgn    = p.Results.IterGaussNewton; 
speak    = p.Results.Verbose; 

C = numel(Nii_x); % Number of channels

msk = cell(1,C); % For saving locations of missing values so that they can be 're-applied' once the algorithm has finished
ll1 = zeros(1,C);
ll2 = 0;
% for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers)  

    % Boundary used to model HR image  
    spm_field('boundary',1);
    pushpull('boundary',1);
    
    % Observed image
    x = get_nii(Nii_x(c));
        
    y = get_nii(Nii_y(c));
    for gnit=1:nitgn % Iterate Gauss-Newton

        % Gradient      
        Ayx  = A(y,dat(c));                
        for n=1:dat(c).N
            % Here we discard missing data, for MRI these are
            % assumed to be zeros and NaNs.
            mskn          = isfinite(x{n}) & x{n} ~= 0;
            Ayx{n}        = Ayx{n} - x{n};
            Ayx{n}(~mskn) = 0;
        end 
        mskn = [];
        rhs  = At(Ayx,dat(c),tau{c})*(1/rho); 
        Ayx  = [];
        rhs  = rhs + spm_field('vel2mom',y,[vx 0 lam(c)^2 0]);

        % Hessian
        H   = get_nii(Nii_H(c));
        lhs = H*sum(tau{c})/rho;
        H   = [];

        % Compute GN step
        y   = y - spm_field(lhs,rhs,[vx 0 lam(c)^2 0 2 2]);
        lhs = [];
        rhs = [];

    end
    msk{c} = isfinite(y) & y ~= 0;       
    
    if strcmpi(modality,'MRI')
        % Ensure non-negativity
        y(y < 0) = 0;
    end 
        
    % Compute log of likelihood part
    ll1(c) = get_ll1(y,x,tau{c},dat(c));
    x      = [];
    
    % Compute log of prior part (part 1)
    G   = lam(c)*imgrad(y,vx);
    ll2 = ll2 + sum(sum(G.^2,4),5);
    G   = [];         
    
    Nii_y(c) = put_nii(Nii_y(c),y);                
    y        = [];
end

% Compute log of prior part (part 2)
ll2 = -sum(sum(sum(sqrt(double(ll2))))); 

if speak >= 2
    % Show initial estimate
    show_progress(method,modality,0,Nii_x,Nii_y,dm); 
end
%==========================================================================
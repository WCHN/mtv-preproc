function [Nii_y,msk] = estimate_initial_y(Nii_x,Nii_y,dat,tau,rho,lam,infnrm,vx,dm,num_workers,p)
% Compute initial estimate of recovered image(s)
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Some parameters from options struct
modality = p.Results.Modality;
method   = p.Results.Method;
nitgn    = p.Results.IterGaussNewton; 

C   = numel(Nii_x);
msk = cell(1,C); % For saving locations of missing values so that they can be 're-applied' once the algorithm has finished
% for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers)  

    spm_field('boundary',1) % Set up boundary conditions that match the gradient operator
            
    % Observed image
    x = get_nii(Nii_x(c));
        
    % Initial guesses for solution    
    if strcmpi(method,'superres')    
        %---------------------------
        % Super-resolution
        %---------------------------
        
        y = get_nii(Nii_y(c));  
        for gnit=1:nitgn % Iterate Gauss-Newton

            % Gradient      
            Ayx  = A(y,dat(c));                
            for n=1:dat(c).N
                % Here we discard missing data, for MRI these are
                % assumed to be zeros and NaNs.
                mskn          = x{n} ~= 0 & isfinite(x{n});
                Ayx{n}        = Ayx{n} - x{n};
                Ayx{n}(~mskn) = 0;
            end 
            mskn = [];
            rhs  = At(Ayx,dat(c),tau{c})*(1/rho); 
            Ayx  = [];
            rhs  = rhs + spm_field('vel2mom',y,[vx 0 lam(c)^2 0]);

            % Hessian
            lhs = infnrm(c)*ones(dm,'single')*sum(tau{c})/rho;

            % Compute GN step
            y   = y - spm_field(lhs,rhs,[vx 0 lam(c)^2 0 2 2]);
            lhs = [];
            rhs = [];

        end
        msk{c} = isfinite(y) & y ~= 0;
    else  
        %---------------------------
        % Denoising
        %---------------------------        
    
        y = spm_field(tau{c}*ones(dm,'single'),tau{c}*x{1},[vx 0 lam(c)^2 0 2 2]); 
        
        if strcmpi(modality,'CT')
            msk{c} = isfinite(x{1}) & x{1} ~= 0 & x{1} ~= min(x{1});
        else
            msk{c} = isfinite(x{1}) & x{1} ~= 0;
        end
    end        
    x = [];
    
    if strcmpi(modality,'MRI')
        % Ensure non-negativity
        y(y < 0) = 0;
    end 
        
    Nii_y(c) = put_nii(Nii_y(c),y);                
    y        = [];
end
%==========================================================================
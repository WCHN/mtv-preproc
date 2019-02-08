function [Nii_y,Nii_u,Nii_w,ll1,ll2]= update_y(Nii_x,Nii_y,Nii_u,Nii_w,dat,tau,rho,lam,infnrm,vx,dm,num_workers,p)
% Compute estimate of recovered image(s)
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Some parameters from options struct
modality = p.Results.Modality;
method   = p.Results.Method;
nitgn    = p.Results.IterGaussNewton; 
speak    = p.Results.Verbose; 
C        = numel(Nii_x);

%------------------------------------------------------------------
% Proximal operator for u
% Here we solve for the MTV term of the objective function using
% vectorial soft thresholding 
%------------------------------------------------------------------

unorm = 0;    
%     for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers) % Loop over channels

    y = get_nii(Nii_y(c));        
    G = lam(c)*imgrad(y,vx);
    y = [];

    w = get_nii(Nii_w(c));        
    u = G + w/rho;
    G = [];
    w = [];

    Nii_u(c) = put_nii(Nii_u(c),u);

    unorm = unorm + sum(sum(u.^2,4),5);
    u     = [];
end % End loop over channels

unorm     = sqrt(unorm);
mtv_scale = max(unorm - 1/rho,0)./(unorm + eps);
clear unorm

ll1 = zeros(1,C);
ll2 = 0;
%     for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers) % Loop over channels

    spm_field('boundary',1) % Set up boundary conditions that match the gradient operator

    u = get_nii(Nii_u(c));   
    w = get_nii(Nii_w(c));   

    %------------------------------------------------------------------
    % Proximal operator for u (continued)
    % Here we multiply each contrast image with the same scaling
    % matrix, this is a key addition of using MTV
    %------------------------------------------------------------------

    u = bsxfun(@times,u,mtv_scale);

    %------------------------------------------------------------------
    % Proximal operator for y        
    %------------------------------------------------------------------

    x = get_nii(Nii_x(c)); % Get observed image

    if strcmpi(method,'superres')              
        %---------------------------
        % Super-resolution
        % Here we want to solve for y using Gauss-Newton (GN)
        % optimisation
        %---------------------------      

        y = get_nii(Nii_y(c)); % Get solution

        for gnit=1:nitgn % Iterate Gauss-Newton

            % Gradient      
            rhs = w/rho - u; 
            rhs = lam(c)*imdiv(rhs,vx);
            Ayx = A(y,dat(c));
            for n=1:dat(c).N
                % Here we discard missing data, for MRI these are
                % assumed to be zeros and NaNs.
                mskn          = isfinite(x{n}) & x{n} ~= 0;
                Ayx{n}        = Ayx{n} - x{n};
                Ayx{n}(~mskn) = 0;
            end                  
            mskn = [];
            rhs  = rhs + At(Ayx,dat(c),tau{c})*(1/rho); 
            Ayx  = [];
            rhs  = rhs + spm_field('vel2mom',y,[vx 0 lam(c)^2 0]);

            % Hessian
            lhs = infnrm(c)*ones(dm,'single')*sum(tau{c})/rho;

            % Compute GN step
            y   = y - spm_field(lhs,rhs,[vx 0 lam(c)^2 0 2 2]);
            lhs = [];
            rhs = [];
        end
    else            
        %---------------------------
        % Denoising
        %---------------------------

        % RHS
        rhs = u - w/rho; 
        rhs = lam(c)*imdiv(rhs,vx);
        rhs = rhs + x{1}*(tau{c}/rho);

        % LHS
        lhs = ones(dm,'single')*tau{c}/rho;

        % Compute new y
        y   = spm_field(lhs,rhs,[vx 0 lam(c)^2 0 2 2]);
        lhs = [];
        rhs = [];
    end                                       

    if strcmpi(modality,'MRI')
        % Ensure non-negativity (ad-hoc)
        y(y < 0) = 0;
    end 

    % Compute likelihood term of posterior
    ll1(c) = get_ll(method,y,x,tau{c},dat(c));
    x      = [];

    %------------------------------------------------------------------
    % Solve for w
    % Here we update the Lagrange variable
    %------------------------------------------------------------------

    G = lam(c)*imgrad(y,vx);

    Nii_y(c) = put_nii(Nii_y(c),y);
    y        = [];

    w = w + rho*(G - u);        

    % Compute prior term of posterior
    ll2 = ll2 + sum(sum(G.^2,4),5);
    G   = [];                

    Nii_u(c) = put_nii(Nii_u(c),u);
    u = [];

    Nii_w(c) = put_nii(Nii_w(c),w);
    w = [];
end % End loop over channels     

if speak >= 2
    % Show MTV scaling
    show_mtv_scale(mtv_scale);
end

% Compute prior part of objective function
ll2 = sum(sum(sum(sqrt(ll2)))); 
%==========================================================================
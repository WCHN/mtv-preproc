function [tau,lam,rho,sched_lam,lam0,Nii_x0] = estimate_model_hyperpars(Nii_x,dec_reg,nit,p)
% Estimate MTV model parameters
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Some parameters from options struct
speak    = p.Results.Verbose; 
rho      = p.Results.ADMMStepSize; 
modality = p.Results.Modality;
method   = p.Results.Method;
if strcmpi(method,'denoise')
    %---------------------------
    % Denoising
    %---------------------------
    scl_lam = p.Results.RegScaleDenoisingMRI;
    lam_ct  = p.Results.RegDenoisingCT; 
elseif strcmpi(method,'superres')
    %---------------------------
    % Super-resolution
    %---------------------------
    scl_lam = p.Results.RegScaleSuperResMRI;
    lam_ct  = p.Results.RegSuperresCT;
end

% Number of channels
C = numel(Nii_x);

% For verbose
Nii_x0 = [];
if speak >= 2
    if strcmpi(modality,'MRI')
        figname = '(SPM) Rice mixture fits MRI';
    elseif strcmpi(modality,'CT')        
        figname = '(SPM) Gaussian mixture fits CT';
    end
    f                = findobj('Type', 'Figure', 'Name', figname);
    if isempty(f), f = figure('Name', figname, 'NumberTitle', 'off'); end
    set(0, 'CurrentFigure', f);              
    
    if speak >= 3
        % So that Verbose = 3 works for superres (because Nii_x are copied, then copies are deleted)
        Nii_x0 = Nii_x; 
    end
end

% Just for making subplots
N0 = 0;
for c=1:C
    N = numel(Nii_x{c});
    for n=1:N
        N0 = N0 + 1;
    end
end
nr = floor(sqrt(N0));
nc = ceil(N0/nr);  

% Make estimates
sd  = cell(1,C);
tau = cell(1,C);
mu  = zeros(1,C);
lam = zeros(1,C);
for c=1:C           
        
    % Estimate image noise and mean brain intensity
    if strcmpi(modality,'MRI')       
        %---------------------------
        % Data is MRI
        %---------------------------
    
        if c > 1, cnt_subplot = cnt_subplot + numel(Nii_x{c - 1});
        else,     cnt_subplot = 0;
        end
        
        [sd{c},mu_brain] = spm_noise_estimate_mod(Nii_x{c},speak >= 2,nr,nc,cnt_subplot); % Noise standard deviation
        
        mu(c)  = mean(mu_brain);        % Mean brain intensity
        lam(c) = scl_lam/double(mu(c)); % This scaling is currently a bit arbitrary, and should be based on empiricism
    elseif strcmpi(modality,'CT')
        %---------------------------
        % Data is CT
        %---------------------------
        
        sd{c} = noise_estimate_ct(Nii_x{c},speak >= 2); % Noise standard deviation
        
        mu(c)  = 0;     % Mean brain intensity not used for CT => intensities follow the Hounsfield scale
        lam(c) = lam_ct;
    end
    
    % Noise precision
    tau{c} = 1./(sd{c}.^2);
end

% For decreasing regularisation with iteration number
lam0      = lam;
sched_lam = get_lam_sched(nit);
if dec_reg
    lam   = sched_lam(1)*lam;
end

if rho == 0
    % Estimate rho (this value seems to lead to reasonably good convergence)
    atau = [];
    for c=1:C
        N = numel(Nii_x{c});
        for n=1:N
            atau = [atau tau{c}(n)];
        end
    end
    rho = sqrt(mean(atau))/mean(lam);    
end

if speak  >= 1
    % Print estimates
    fprintf('Estimated parameters are:\n');
    for c=1:C        
        N = numel(Nii_x{c});
        for n=1:N
            fprintf('c=%i | n=%i | sd=%f, mu=%f | tau=%f, lam=%f, rho=%f\n', c, n, sd{c}(n), mu(c), tau{c}(n), lam(c), rho);
        end
    end
    fprintf('\n');
end
%==========================================================================

%==========================================================================
function sched = get_lam_sched(nit)
def        = spm_shoot_defaults;
sched      = def.sched(1:14);
sched(end) = 1;
% sched = fliplr(1:2:20);
%==========================================================================
function [tau,lam,sched] = estimate_model_hyperpars(Nii_x,dec_reg,vx,modality,p)
% Estimate MTV model parameters
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Some parameters from options struct
speak    = p.Results.Verbose; 
rho      = p.Results.ADMMStepSize; 
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
sd          = cell(1,C);
tau         = cell(1,C);
mu          = zeros(1,C);
lam         = zeros(1,C);
cnt_subplot = 1;
if speak >= 2, clear_figure(modality); end
for c=1:C           
        
    % Estimate image noise and mean brain intensity
    if strcmpi(modality,'MRI')       
        %---------------------------
        % Data is MRI
        %---------------------------
        
        if speak >= 2, set_figure(modality); end
        
        % Compute noise standard deviation            
        try
            [sd{c},mu_brain,info] = spm_noise_estimate(Nii_x{c},3); 
        catch
            [sd{c},mu_brain,info] = my_spm_noise_estimate(Nii_x{c},3);
        end
        
        if speak >= 2
            % Show fit   
            for i=1:numel(info)
                subplot(nr,nc,cnt_subplot);                        
                plot(info(i).x(:),info(i).p,'--',info(i).x(:),info(i).h/sum(info(i).h)/info(i).md,'b.',info(i).x(:),info(i).sp,'r'); drawnow
                cnt_subplot = cnt_subplot + 1;
            end
        end
        
        mu(c)  = mean(mu_brain);        % Mean brain intensity
        lam(c) = scl_lam/double(mu(c)); % This scaling is currently a bit arbitrary, and should be based on empiricism
    elseif strcmpi(modality,'CT')
        %---------------------------
        % Data is CT
        %---------------------------
        
        if speak >= 2, set_figure(modality); end
        
        [sd{c},info] = noise_estimate_ct(Nii_x{c}); % Noise standard deviation
           
        if speak >= 2
            % Show fit
            subplot(121)
            plot(info.x(:),info.p,'--',info.x(:),info.h/sum(info.h)/info.md,'b.',info.x(:),info.sp,'r'); 
            subplot(122)
            hist(info.f(:),info.x')
            drawnow
        end    
    
        mu(c)  = 0;     % Mean brain intensity not used for CT => intensities follow the Hounsfield scale
        lam(c) = lam_ct;
    end

    % Noise precision
    tau{c} = 1./(sd{c}.^2);
end
clear info

% Get all stds
all_sd = [];
for c=1:C        
    N = numel(Nii_x{c});
    for n=1:N
        all_sd = [all_sd sd{c}(n)];
    end
end
    
% Incorporate template voxel size into regularisation
lam = (prod(vx))^(1/2)*lam;

% % An attempt to combat the bias-variance trade-off -> scale lambda by the 
% % number of observations of each channel
% for c=1:C
%     N      = numel(tau{c});
%     lam(c) = sqrt(N)*lam(c);
% end

% For decreasing regularisation with iteration number
sched = get_lam_sched(mu,tau,scl_lam,all_sd,dec_reg);

if rho == 0
    % Estimate rho (this value seems to lead to reasonably good convergence)
    rho = estimate_rho(tau,lam);       
end
    
if speak  >= 1
    % Print estimates
    fprintf('Estimated parameters are:\n');
    for c=1:C        
        N = numel(Nii_x{c});
        for n=1:N
            fprintf('c=%1i | n=%1i | sd=%10.4f, mu=%10.4f | tau=%10.4e, lam=%10.4e, rho=%10.4g\n', c, n, sd{c}(n), mu(c), tau{c}(n), lam(c), rho);
        end
    end
    fprintf('\n');
end
%==========================================================================

%==========================================================================
function sched = get_lam_sched(mu,tau,scl_lam,asd,dec_reg)
C         = numel(tau);
sched     = struct;
sched.it  = 1;
sched.cnt = 1;    
if dec_reg
    vals      = exp(linspace(log(2^5),0,6));
%     vals      = [16 8 6 4 2 1];
%     vals      = linspace(50,0,10);
    sched.scl = zeros(numel(vals),C);
    for c=1:C
         sched.scl(:,c)  = vals;
    end    
    sched.nxt = 1;
else
    sched.scl = ones(1,C);
    sched.nxt = Inf;
end
%==========================================================================

%==========================================================================
function f = set_figure(modality)
if strcmpi(modality,'MRI')
    figname = '(SPM) Rice mixture fits MRI';
elseif strcmpi(modality,'CT')        
    figname = '(SPM) Gaussian mixture fits CT';
end
f                = findobj('Type', 'Figure', 'Name', figname);
if isempty(f), f = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', f);    
%==========================================================================    

%==========================================================================
function clear_figure(modality)
fig = set_figure(modality);
clf(fig)
%==========================================================================    
function Nio = MTVprocess3D(varargin)
% Multi-channel total variation (MTV) denoising or super-resolution of MR images. 
% Requires that the SPM software is on the MATLAB path. SPM is available from:
% https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
%
% FORMAT Nio = MTVprocess3D(...)
%
% KEYWORD
% -------
%
% InputImages          - Either image filenames in a cell array, or images in a nifti object. If empty, uses spm_select ['']
% IterMax              - Maximum number of iteration [30]
% Tolerance            - Convergence threshold [1e-3]
% Regularisation_scale - Scaling of regularisation, increase this value for stronger denoising [15]
% WorkersParfor        - Maximum number of parfor workers [Inf]
% Temporary_directory  - Directory for temporary files ['./tmp']
% Output_directory     - Directory for denoised images ['./out']
% Method               - Does either denoising ('denoise') or super-resolution ('superres') ['denoise']
% Verbose              - Verbosity level:  0  = quiet
%                                         [1] = write  (log likelihood, parameter estimates)
%                                          2  = draw   (log likelihood, rice fit, noisy+cleaned)
%                                          3  = result (show noisy and denoised image(s) in spm_check_registration)
% CleanUp              - Delete temporary files [true] 
% 
% OUTPUT
% ------
% 
% Nio - nifti object containing denoised images
% 
%__________________________________________________________________________
% The general principles are described in the following paper:
%
%     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
%     MRI Super-Resolution Using Multi-channel Total Variation.
%     In Annual Conference on Medical Image Understanding and Analysis
%     2018 Jul 9 (pp. 217-228). Springer, Cham.
%
% OBS: The code uses MATLAB's parfor to parallelise and speed up certain
% processing. The code should be memory efficient, still, running parfor
% can lead to the use of more RAM than what is available. To decrease the
% number of parfor workers, use the WorkersParfor option described below.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
% Parse options
%--------------------------------------------------------------------------

p              = inputParser;
p.FunctionName = 'MTVprocess3D';
p.addParameter('InputImages', '', @(in) (ischar(in) || isa(in,'nifti')));
p.addParameter('IterMax', 30, @isnumeric);
p.addParameter('Tolerance', 1e-3, @isnumeric);
p.addParameter('Regularisation_scale', 15, @isnumeric);
p.addParameter('WorkersParfor', Inf, @(in) (isnumeric(in) && in >= 0));
p.addParameter('Temporary_directory', 'tmp', @ischar);
p.addParameter('Output_directory', 'out', @ischar);
p.addParameter('Method', 'denoise', @(in) (ischar(in) && (strcmpi(in,'denoise') || strcmpi(in,'superres'))));
p.addParameter('Verbose', 1, @(in) (isnumeric(in) && in >= 0 && in <= 3));
p.addParameter('CleanUp', true, @islogical);
p.parse(varargin{:});
nii_x       = p.Results.InputImages;
nit         = p.Results.IterMax;
tol         = p.Results.Tolerance;
scl_lam     = p.Results.Regularisation_scale;
num_workers = p.Results.WorkersParfor;
dir_tmp     = p.Results.Temporary_directory;
dir_out     = p.Results.Output_directory;
method      = p.Results.Method;
speak       = p.Results.Verbose; 
do_clean    = p.Results.CleanUp; 

if strcmpi(method,'superres'), error('Super-resolution not yet added to this code. But coming soon!'); end

% Make some directories
if  exist(dir_tmp,'dir'), rmdir(dir_tmp,'s'); end; mkdir(dir_tmp); 
if ~exist(dir_out,'dir'), mkdir(dir_out);  end

% Add to MATLAB path
addpath(fullfile(fileparts(mfilename('fullpath')),'code'));
  
%--------------------------------------------------------------------------
% Get image data
%--------------------------------------------------------------------------

if isempty(nii_x)
    nii_x = nifti(spm_select(Inf,'nifti','Select image'));    
else
    if ~isa(nii_x,'nifti'), nii_x = nifti(nii_x); end
end
C = numel(nii_x); % Number of channels

% Sanity check input
for c=1:C    
    dm = nii_x(c).dat.dim;
    
    if c > 1 && (~isequal(dm,odm) || ~isequal(dm,odm))
        error('Images are not all the same size!')
    end
    odm = dm;
end
mat = nii_x(1).mat;
vx  = sqrt(sum(nii_x(c).mat(1:3,1:3).^2));

%--------------------------------------------------------------------------
% Estimate model hyper-parameters
%--------------------------------------------------------------------------

if speak >= 2
    figname          = '(SPM) Rice mixture fits';
    f                = findobj('Type', 'Figure', 'Name', figname);
    if isempty(f), f = figure('Name', figname, 'NumberTitle', 'off'); end
    set(0, 'CurrentFigure', f);  
    
    nr = floor(sqrt(C));
    nc = ceil(C/nr);  
end

sd  = zeros(1,C);
tau = zeros(1,C);
mu  = zeros(1,C);
lam = zeros(1,C);
for c=1:C           
    if speak >= 2, subplot(nr,nc,c); end
    
    % Estimate image noise and mean brain intensity
    [sd(c),mu_brain] = spm_noise_estimate_mod(nii_x(c),speak >= 2);
    
    tau(c) = 1/(sd(c).^2);          % Noise precision
    mu(c)  = mu_brain;              % Mean brain intensity
    lam(c) = scl_lam/double(mu(c)); % This scaling is currently a bit arbitrary, and should be based on empiricism
end

% Estimate rho
rho = sqrt(mean(tau))/mean(lam); % This value of rho seems to lead to reasonably good convergence

if speak  >= 1
    % Print estimates
    fprintf('Estimated parameters are:\n');
    for c=1:C
        fprintf('c=%i -> sd=%f, mu=%f -> tau=%f, lam=%f, rho=%f\n', c, sd(c), mu(c), tau(c), lam(c), rho);
    end
    fprintf('\n');
end

% Set up boundary condidtions that match the gradient operator
spm_field('boundary',1)

%--------------------------------------------------------------------------
% Initialise NIfTIs
%--------------------------------------------------------------------------

nii_y = nifti;
nii_u = nifti;
nii_w = nifti;
for c=1:C
    fname_y = fullfile(dir_tmp,['y' num2str(c) '.nii']);
    fname_u = fullfile(dir_tmp,['u' num2str(c) '.nii']); 
    fname_w = fullfile(dir_tmp,['w' num2str(c) '.nii']);
    
    spm_misc('create_nii',fname_y,zeros(dm,'single'),mat,[spm_type('float32') spm_platform('bigend')],'y');
    spm_misc('create_nii',fname_u,zeros([dm 3],'single'),mat,[spm_type('float32') spm_platform('bigend')],'u');
    spm_misc('create_nii',fname_w,zeros([dm 3],'single'),mat,[spm_type('float32') spm_platform('bigend')],'w');
    
    nii_y(c) = nifti(fname_y);
    nii_u(c) = nifti(fname_u);
    nii_w(c) = nifti(fname_w);
end

%--------------------------------------------------------------------------
% Initialise variables
%--------------------------------------------------------------------------

% Manage parfor
if num_workers == Inf, num_workers = spm_misc('nbr_parfor_workers'); end
spm_misc('manage_parpool',min(C,num_workers));

ll1 = zeros(1,C);
ll2 = 0;
parfor c=1:C
    
    % Noisy image
    x = get_nii(nii_x,c);    
    
    % Denoised image
    y = spm_field(tau(c)*ones(dm,'single'),tau(c)*x,[vx  0 lam(c)^2 0  2 2]); 
    put_nii(nii_y,c,y);
        
    % Objective
    ll1(c) = -(tau(c)/2)*sum(sum(sum((y - x).^2)));
    x      = [];         
    
    G = imgrad(y,lam(c),dm,vx);
    y = [];
    
    ll2 = ll2 + sum(G.^2,4);
    G   = []; 
         
    % Proximal variables
    u = zeros([dm 3],'single');
    put_nii(nii_u,c,u);
    u = [];
    
    w = zeros([dm 3],'single');    
    put_nii(nii_w,c,w);
    w = [];
end

% Objective
ll2 = -sum(sum(sum(sqrt(ll2))));
ll  = -sum(ll1) + ll2;

%--------------------------------------------------------------------------
% Start denoising
%--------------------------------------------------------------------------

if speak >= 1, fprintf('Started method: %s --- running (max) %d iterations\n', method, nit); end

if speak >= 1, tic; end
for it=1:nit
        
    %------------------------------------------------------------------
    % Proximal operator for u
    %------------------------------------------------------------------

    unorm = 0;    
    parfor c=1:C  % Loop over channels
        y = get_nii(nii_y,c);        
        G = imgrad(y,lam(c),dm,vx);
        y = [];
             
        w = get_nii(nii_w,c);        
        u = G + w/rho;
        G = [];
        w = [];
        
        put_nii(nii_u,c,u);
        
        unorm = unorm + sum(u.^2,4);
        u     = [];
    end
    unorm = sqrt(unorm);
    scale = max(unorm - 1/rho,0)./(unorm + eps);
    clear unorm
    
    ll1 = zeros(1,C);
    ll2 = 0;
    parfor c=1:C  % Loop over channels
                
        u = get_nii(nii_u,c);   
        w = get_nii(nii_w,c);   
        
        %------------------------------------------------------------------
        % Proximal operator for u (continued)
        %------------------------------------------------------------------
        
        u = bsxfun(@times,u,scale);
            
        %------------------------------------------------------------------
        % Proximal operator for y
        %------------------------------------------------------------------
                      
        g = u - w/rho; 
        g = imdiv(g,lam(c),dm,vx);
                
        x = get_nii(nii_x,c);    
        g = g + x*(tau(c)/rho);
        
        y = spm_field(ones(dm,'single')*tau(c)/rho,g,[vx  0 lam(c)^2 0  2 2]);
        g = [];

        % Objective function
        ll1(c) = -(tau(c)/2)*sum(sum(sum((y - x).^2)));
        x      = [];
        
        %------------------------------------------------------------------
        % Solve for w
        %------------------------------------------------------------------
        
        G = imgrad(y,lam(c),dm,vx);
        
        put_nii(nii_y,c,y);
        y = [];
        
        w = w - rho*(u - G);        
    
        % Objective function    
        ll2 = ll2 + sum(G.^2,4);      
        G   = [];                
        
        put_nii(nii_u,c,u);
        u = [];
        
        put_nii(nii_w,c,w);
        w = [];
    end        
    clear scale
    
    % Objective function    
    ll2  = -sum(sum(sum(sqrt(ll2))));    
    ll   = [ll, sum(ll1) + ll2];    
    gain = abs((ll(end - 1)*(1 + 10*eps) - ll(end))/ll(end));
    
    % Some (potential) verbose            
    if speak >= 1, fprintf('%d %g %g %g %g %g\n', it, sum(ll1), ll2, sum(ll1) + ll2,gain,tol); end
    if speak >= 2, show_progress(ll,nii_x,nii_y,dm,nr,nc); end
    
    if gain < tol && it > 10
        % Finished        
        break;
    end
end
clear u w G x D

if speak >= 1, toc; end

%--------------------------------------------------------------------------
% Write results
%--------------------------------------------------------------------------

Nio = nifti;
for c=1:C
    Nio(c)           = nii_x(c);
    [~,nam,ext]      = fileparts(Nio(c).dat.fname);
    Nio(c).dat.fname = fullfile(dir_out,['den_' nam ext]);
    create(Nio(c));
    
    y                 = get_nii(nii_y,c);        
    Nio(c).dat(:,:,:) = y;
end

%--------------------------------------------------------------------------
% Show original and denoised
%--------------------------------------------------------------------------

if speak >= 3
    fnames = cell(1,2*C);
    cnt    = 1;
    for c=1:2:2*C    
        fnames{c}     = nii_x(cnt).dat.fname;    
        fnames{c + 1} = Nio(cnt).dat.fname;
        cnt           = cnt + 1;
    end

    spm_check_registration(char(fnames))
end

if do_clean
    % Clean-up
    rmdir(dir_tmp,'s');
end
%==========================================================================
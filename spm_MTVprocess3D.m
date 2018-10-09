function Nio = spm_MTVprocess3D(varargin)
% Multi-channel total variation (MTV) denoising or super-resolution of MR images. 
% Requires that the SPM software is on the MATLAB path. SPM is available from:
% https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
%
% FORMAT Nio = MTVprocess3D(...)
%
% KEYWORD
% -------
%
% InputImages          - Either image filenames in a cell array, or images 
%                        in a nifti object. If empty, uses spm_select ['']
% IterMax              - Maximum number of iteration [30]
% Tolerance            - Convergence threshold [1e-3]
% Regularisation_scale - Scaling of regularisation, increase this value for 
%                        stronger denoising [15]
% WorkersParfor        - Maximum number of parfor workers [Inf]
% TemporaryDirectory  - Directory for temporary files ['./tmp']
% OutputDirectory     - Directory for denoised images ['./out']
% Method               - Does either denoising ('denoise') or 
%                        super-resolution ('superres') ['denoise']
% Verbose              - Verbosity level:  0  = quiet
%                                         [1] = write  (log likelihood, parameter estimates)
%                                          2  = draw   (log likelihood, rice fit, noisy+cleaned)
%                                          3  = result (show noisy and denoised image(s) in spm_check_registration)
% CleanUp              - Delete temporary files [true] 
% VoxelSize            - Voxel size of super-resolved image [1 1 1]
% IterMaxCG            - Maximum number of iterations for conjugate gradient 
%                        solver used for super-resolution [20]
% ToleranceCG          - Convergence threshold for conjugate gradient 
%                        solver used for super-resolution [1e-3]
% CoRegister           - For super-resolution, co-register input images [true] 
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

% First check that SPM is on the MATLAB path
if ~(exist('spm','file') == 2), error('SPM is not on the MATLAB path!'); end

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
p.addParameter('TemporaryDirectory', 'tmp', @ischar);
p.addParameter('OutputDirectory', 'out', @ischar);
p.addParameter('Method', 'denoise', @(in) (ischar(in) && (strcmpi(in,'denoise') || strcmpi(in,'superres'))));
p.addParameter('Verbose', 1, @(in) (isnumeric(in) && in >= 0 && in <= 3));
p.addParameter('CleanUp', true, @islogical);
p.addParameter('VoxelSize', [1 1 1], @(in) (isnumeric(in) && (numel(in) == 1 || numel(in) == 3)) && ~any(in <= 0));
p.addParameter('IterMaxCG', 4, @isnumeric);
p.addParameter('ToleranceCG', 1e-3, @isnumeric);
p.addParameter('CoRegister', true, @islogical);
p.parse(varargin{:});
Nii_x       = p.Results.InputImages;
nit         = p.Results.IterMax;
tol         = p.Results.Tolerance;
scl_lam     = p.Results.Regularisation_scale;
num_workers = p.Results.WorkersParfor;
dir_tmp     = p.Results.TemporaryDirectory;
dir_out     = p.Results.OutputDirectory;
method      = p.Results.Method;
speak       = p.Results.Verbose; 
do_clean    = p.Results.CleanUp; 
vx_sr       = p.Results.VoxelSize; 
nit_cg      = p.Results.IterMaxCG; 
tol_cg      = p.Results.ToleranceCG; 
do_coreg    = p.Results.CoRegister; 

% Make some directories
if  exist(dir_tmp,'dir'), rmdir(dir_tmp,'s'); end; mkdir(dir_tmp); 
if ~exist(dir_out,'dir'), mkdir(dir_out);  end
  
% Set up boundary conditions that match the gradient operator
spm_field('boundary',1)

%--------------------------------------------------------------------------
% Get image data
%--------------------------------------------------------------------------

if isempty(Nii_x)
    Nii_x = nifti(spm_select(Inf,'nifti','Select image'));    
else
    if ~isa(Nii_x,'nifti'), Nii_x = nifti(Nii_x); end
end
Nii_x0 = Nii_x;        % So that Verbose = 3 works for superres
C      = numel(Nii_x); % Number of channels

% Sanity check input
for c=1:C    
    dm = Nii_x(c).dat.dim;
    
    if strcmpi(method,'denoise') && c > 1 && (~isequal(dm,odm) || ~isequal(dm,odm))
        error('Images are not all the same size!')
    end
    odm = dm;
end

% Get voxel size, orientation matrix and image dimensions
if strcmpi(method,'denoise')
    mat = Nii_x(1).mat;
    vx  = sqrt(sum(mat(1:3,1:3).^2));    
elseif strcmpi(method,'superres')
    % For super-resolution, calculate orientation matrix and dimensions 
    % from maximum bounding-box
    vx       = vx_sr;
    [mat,dm] = max_bb_orient(Nii_x,vx);
end

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
    [sd(c),mu_brain] = spm_noise_estimate_mod(Nii_x(c),speak >= 2);
    
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

%--------------------------------------------------------------------------
% Initialise NIfTIs
%--------------------------------------------------------------------------

Nii_y = nifti;
Nii_u = nifti;
Nii_w = nifti;
for c=1:C
    fname_y = fullfile(dir_tmp,['y' num2str(c) '.nii']);
    fname_u = fullfile(dir_tmp,['u' num2str(c) '.nii']); 
    fname_w = fullfile(dir_tmp,['w' num2str(c) '.nii']);
    
    create_nii(fname_y,zeros(dm,    'single'),mat,[spm_type('float32') spm_platform('bigend')],'y');
    create_nii(fname_u,zeros([dm 3],'single'),mat,[spm_type('float32') spm_platform('bigend')],'u');
    create_nii(fname_w,zeros([dm 3],'single'),mat,[spm_type('float32') spm_platform('bigend')],'w');
    
    Nii_y(c) = nifti(fname_y);
    Nii_u(c) = nifti(fname_u);
    Nii_w(c) = nifti(fname_w);
end

%--------------------------------------------------------------------------
% Initialise variables
%--------------------------------------------------------------------------

if strcmpi(method,'superres')    
    if do_coreg
        % Co-register input images
        Nii_x = coreg(Nii_x,dir_tmp);
    end
    
    % Initialise dat struct with projection matrices, etc.
    dat = init_dat(Nii_x,mat,dm);  
    
    % Define Laplace prior (in Fourier space)
    L = laplace_prior(dm,vx);
else 
    % For denoising these variables need to be defined, otherwise parfor
    % will complain :(
    dat = zeros(1,C);
    L   = [];
end

% Manage parfor
if num_workers == Inf, num_workers = nbr_parfor_workers; end
manage_parpool(min(C,num_workers));

parfor c=1:C    
    
    % Noisy image
    x = get_nii(Nii_x(c));    
    
    % Initial guess for solution    
    if strcmpi(method,'superres')    
        % Super-resolved image (using bspline interpolation)
        y = get_y_superres(Nii_x(c),dat(c),dm,mat);
    else  
        % Denoised image
        y = spm_field(tau(c)*ones(dm,'single'),tau(c)*x,[vx  0 lam(c)^2 0  2 2]); 
    end
    put_nii(Nii_y(c),y);                
    y = [];
    
    % Proximal variables
    u = zeros([dm 3],'single');
    put_nii(Nii_u(c),u);
    u = [];
    
    w = zeros([dm 3],'single');    
    put_nii(Nii_w(c),w);
    w = [];
end

%--------------------------------------------------------------------------
% Start denoising
%--------------------------------------------------------------------------

if speak >= 1, fprintf('Started method: %s --- running (max) %d iterations\n', method, nit); end

if speak >= 1, tic; end
for it=1:nit
        
    %------------------------------------------------------------------
    % Proximal operator for u
    % Here we solve for the MTV term of the objective function using the
    % vectorial soft thresholding 
    %------------------------------------------------------------------

    unorm = 0;    
    parfor c=1:C  % Loop over channels
        
        y = get_nii(Nii_y(c));        
        G = lam(c)*imgrad(y,vx);
        y = [];
             
        w = get_nii(Nii_w(c));        
        u = G + w/rho;
        G = [];
        w = [];
        
        put_nii(Nii_u(c),u);
        
        unorm = unorm + sum(u.^2,4);
        u     = [];
    end
    unorm = sqrt(unorm);
    scale = max(unorm - 1/rho,0)./(unorm + eps);
    clear unorm
        
    ll1 = zeros(1,C);
    ll2 = 0;
    parfor c=1:C  % Loop over channels
                
        u = get_nii(Nii_u(c));   
        w = get_nii(Nii_w(c));   
        
        %------------------------------------------------------------------
        % Proximal operator for u (continued)
        % Here we multiply each contrast image with the same scaling
        % matrix, this is a key addition of using MTV
        %------------------------------------------------------------------
        
        u = bsxfun(@times,u,scale);
            
        %------------------------------------------------------------------
        % Proximal operator for y
        % Here we want to solve the linear system: lhs\rhs
        %------------------------------------------------------------------
                  
        rhs = u - w/rho; 
        rhs = lam(c)*imdiv(rhs(:,:,:,1),rhs(:,:,:,2),rhs(:,:,:,3),vx);
               
        x = get_nii(Nii_x(c));    
        
        if strcmpi(method,'superres')  
            % Super-resolution
            rhs = rhs + At({x},tau(c),dat(c))*(1/rho);            
            lhs = @(y) AtA(y,@(y) L(y),tau(c)/rho,lam(c)^2,dat(c));
            
            % Solve using conjugate gradient
            y = cg_im_solver(lhs,rhs,get_nii(Nii_y(c)),nit_cg,tol_cg);
            
            % Ensure non-negative
            y(y < 0) = 0;
        else            
            % Denoising
            rhs = rhs + x*(tau(c)/rho);
            lhs = ones(dm,'single')*tau(c)/rho;
            
            % Solve using full multi-grid
            y = spm_field(lhs,rhs,[vx  0 lam(c)^2 0  2 2]);
        end                               
        lhs = [];
        rhs = [];

        % Compute likelihood term of posterior
        ll1(c) = get_ll1(method,y,x,tau(c),dat(c));
        x      = [];
        
        %------------------------------------------------------------------
        % Solve for w
        % Here we update the Lagrange variable
        %------------------------------------------------------------------
        
        G = lam(c)*imgrad(y,vx);
        
        put_nii(Nii_y(c),y);
        y = [];
        
        w = w - rho*(u - G);        
    
        % Compute prior term of posterior
        ll2 = ll2 + sum(G.^2,4);      
        G   = [];                
        
        put_nii(Nii_u(c),u);
        u = [];
        
        put_nii(Nii_w(c),w);
        w = [];
    end        
    clear scale
    
    % Compute log-posterior (objective value)
    ll2  = -sum(sum(sum(sqrt(ll2))));    
    ll   = [ll, sum(ll1) + ll2];    
    gain = abs((ll(end - 1)*(1 + 10*eps) - ll(end))/ll(end));
    
    % Some (potential) verbose            
    if speak >= 1, fprintf('%d %g %g %g %g %g\n', it, sum(ll1), ll2, sum(ll1) + ll2,gain,tol); end
    if speak >= 2, show_progress(ll,Nii_x,Nii_y,dm,nr,nc); end
    if speak >= 2, show_progress(method,ll,Nii_x,Nii_y,dm,nr,nc); end
    
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

if strcmpi(method,'superres'), prefix = 'sr';
else,                          prefix = 'den';
end
   
Nii = nifti;
for c=1:C
    VO           = spm_vol(Nii_x(c).dat.fname);
    [~,name,ext] = fileparts(VO.fname);
    VO.fname     = fullfile(dir_out,[prefix '_' name ext]);
    VO.dim(1:3)  = dm(1:3);    
    VO.mat       = mat;
    VO           = spm_create_vol(VO);
        
    Nii(c)            = nifti(VO.fname);
    y                 = get_nii(Nii_y(c)); 
    Nii(c).dat(:,:,:) = y;    
end

%--------------------------------------------------------------------------
% Show input and solved
%--------------------------------------------------------------------------

if speak >= 3
    fnames = cell(1,2*C);
    cnt    = 1;
    for c=1:2:2*C    
        fnames{c}     = Nii_x0(cnt).dat.fname;    
        fnames{c + 1} = Nii(cnt).dat.fname;
        cnt           = cnt + 1;
    end

    spm_check_registration(char(fnames))
end

if do_clean
    % Clean-up
    rmdir(dir_tmp,'s');
end
%==========================================================================
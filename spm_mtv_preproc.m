function Nii = spm_mtv_preproc(varargin)
% Multi-channel total variation (MTV) denoising or super-resolution of 
% MR images. 
%
% Requires that the SPM software is on the MATLAB path.
% SPM is available from: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
%
% FORMAT Nio = spm_mtv_preproc(...)
%
% KEYWORD
% -------
%
% InputImages            - Either image filenames in a cell array, or images 
%                          in a nifti object. If empty, uses spm_select ['']
% IterMax                - Maximum number of iteration [40]
% ADMMStepSize           - The infamous ADMM step size, set to zero for an 
%                          educated guess [0.1]
% Tolerance              - Convergence threshold, set to zero to run until 
%                          IterMax [0]
% RegularisationScaleMRI - Scaling of regularisation, increase this value for 
%                          stronger denoising [20]
% WorkersParfor          - Maximum number of parfor workers [Inf]
% TemporaryDirectory     - Directory for temporary files ['./tmp']
% OutputDirectory        - Directory for denoised images ['./out']
% Method                 - Does either denoising ('denoise') or 
%                          super-resolution ('superres') ['denoise']
% Verbose                - Verbosity level: 
%                          *  0  = quiet
%                          * [1] = write  (log likelihood, parameter estimates)
%                          *  2  = draw   (log likelihood, rice fit, noisy+cleaned)
%                          *  3  = result (show noisy and denoised image(s) in spm_check_registration)
% CleanUp                - Delete temporary files [true] 
% VoxelSize              - Voxel size of super-resolved image [1 1 1]
% IterMaxCG              - Maximum number of iterations for conjugate gradient 
%                          solver used for super-resolution [12]
% ToleranceCG            - Convergence threshold for conjugate gradient 
%                          solver used for super-resolution [1e-3]
% CoRegister             - For super-resolution, co-register input images [true] 
% Modality               - Either MRI (denoise and super-resolution) or CT 
%                          (denoise) ['MRI']
% RegularisationCT       - Regularisation used for CT denoising [0.04]
% ReadWrite              - Keep variables in workspace (requires more RAM,
%                          but faster), or read/write from disk (requires 
%                          less RAM, but slower) [false] 
% SuperResWithFMG        - Use either spm_field (true) or conjugate
%                          gradient (false) for super-resolution [true]
% ZeroMissingValues      - Set NaNs and zero values to zero after algorithm 
%                          has finished [C=1:true, C>1:false]
%
% OUTPUT
% ------
% 
% Nii - nifti object containing denoised/super-resolved images
% 
%__________________________________________________________________________
%
% Example 1: Super-resolve a set MRIs of one subject
%
% Generate thick-sliced from IXI references by running the script:
% >> GenerateTestData % Down-sampling factor set by DownSampling parameter
%
% Read thick-sliced IXI MRIs
% >> dir_data = './data';
% >> Nii      = nifti(spm_select('FPList',dir_data,'^ds_.*\.nii$'));
%
% Super-resolve the MRIs
% >> spm_mtv_preproc('InputImages',Nii,'Method','superres','Verbose',2);
%
% Compare super-resolved with known ground-truth
% >> files_sr  = spm_select('FPList','./out', '^sr_.*\.nii$');
% >> files_ref = spm_select('FPList',dir_data,'^IXI.*\.nii$');
% >> spm_check_registration(char({files_sr,files_ref}));
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
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% First check that all is okay with SPM
spm_check_path('pull');

%--------------------------------------------------------------------------
% Parse input
%--------------------------------------------------------------------------

p              = inputParser;
p.FunctionName = 'spm_mtv_preproc';
p.addParameter('InputImages', '', @(in) (ischar(in) || isa(in,'nifti')));
p.addParameter('IterMax', 40, @(in) (isnumeric(in) && in > 0));
p.addParameter('ADMMStepSize', 0.1, @(in) (isnumeric(in) && in >= 0));
p.addParameter('Tolerance', 0, @(in) (isnumeric(in) && in >= 0));
p.addParameter('RegularisationScaleMRI', 20, @(in) (isnumeric(in) && in > 0));
p.addParameter('WorkersParfor', Inf, @(in) (isnumeric(in) && in >= 0));
p.addParameter('TemporaryDirectory', 'tmp', @ischar);
p.addParameter('OutputDirectory', 'out', @ischar);
p.addParameter('Method', 'denoise', @(in) (ischar(in) && (strcmpi(in,'denoise') || strcmpi(in,'superres'))));
p.addParameter('Verbose', 1, @(in) (isnumeric(in) && in >= 0 && in <= 3));
p.addParameter('CleanUp', true, @islogical);
p.addParameter('VoxelSize', [1 1 1], @(in) (isnumeric(in) && (numel(in) == 1 || numel(in) == 3)) && ~any(in <= 0));
p.addParameter('IterMaxCG', 12, @(in) (isnumeric(in) && in > 0));
p.addParameter('ToleranceCG', 1e-4, @(in) (isnumeric(in) && in >= 0));
p.addParameter('CoRegister', true, @islogical);
p.addParameter('Modality', 'MRI', @(in) (ischar(in) && (strcmpi(in,'MRI') || strcmpi(in,'CT'))));
p.addParameter('RegularisationCT', 0.04, @(in) (isnumeric(in) && in > 0));
p.addParameter('ReadWrite', false, @islogical);
p.addParameter('SuperResWithFMG', true, @islogical);
p.addParameter('ZeroMissingValues', [], @(in) (islogical(in) || isnumeric(in)));
p.parse(varargin{:});
Nii_x        = p.Results.InputImages;
nit          = p.Results.IterMax;
tol          = p.Results.Tolerance;
scl_lam      = p.Results.RegularisationScaleMRI;
num_workers  = p.Results.WorkersParfor;
dir_tmp      = p.Results.TemporaryDirectory;
dir_out      = p.Results.OutputDirectory;
method       = p.Results.Method;
speak        = p.Results.Verbose; 
do_clean     = p.Results.CleanUp; 
vx_sr        = p.Results.VoxelSize; 
nit_cg       = p.Results.IterMaxCG; 
tol_cg       = p.Results.ToleranceCG; 
coreg        = p.Results.CoRegister; 
modality     = p.Results.Modality; 
lam_ct       = p.Results.RegularisationCT; 
do_readwrite = p.Results.ReadWrite; 
superes_fmg  = p.Results.SuperResWithFMG; 
rho          = p.Results.ADMMStepSize; 
zeroMissing  = p.Results.ZeroMissingValues; 

if strcmpi(method,'superres') && strcmpi(modality,'CT')
    error('Super-resolution not yet supported for CT data!');
end
  
if numel(vx_sr) == 1, vx_sr = vx_sr*ones(1,3); end

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

if isempty(zeroMissing)    
    % Missing values (NaNs and zeros) will be...
    if C == 1
        % ...set to zero after algorithm finishes, if only ONE channel
        zeroMissing = true;
    else
        % ...filled in by the algorithm, if MORE than one channel
        zeroMissing = false;
    end
end

% Make some directories
if  exist(dir_tmp,'dir') == 7,  rmdir(dir_tmp,'s'); end
if  do_readwrite || (coreg && C > 1), mkdir(dir_tmp); end
if ~(exist(dir_out,'dir') == 7),  mkdir(dir_out);  end

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
    %---------------------------
    % Denoising
    %---------------------------
    
    mat = Nii_x(1).mat;
    vx  = sqrt(sum(mat(1:3,1:3).^2));   
elseif strcmpi(method,'superres')
    %---------------------------
    % Super-resolution
    %---------------------------
            
    % For super-resolution, calculate orientation matrix and dimensions 
    % from maximum bounding-box
    vx          = vx_sr;
    [mat,dm]    = max_bb_orient(Nii_x,vx);
    zeroMissing = false;
end

%--------------------------------------------------------------------------
% Estimate model hyper-parameters
%--------------------------------------------------------------------------

if speak >= 2
    if strcmpi(modality,'MRI')
        figname          = '(SPM) Rice mixture fits';
        f                = findobj('Type', 'Figure', 'Name', figname);
        if isempty(f), f = figure('Name', figname, 'NumberTitle', 'off'); end
        set(0, 'CurrentFigure', f);  
    end
    
    nr = floor(sqrt(C));
    nc = ceil(C/nr);  
end

sd  = zeros(1,C);
tau = zeros(1,C);
mu  = zeros(1,C);
lam = zeros(1,C);
for c=1:C           
    if speak >= 2 && strcmpi(modality,'MRI'), subplot(nr,nc,c); end
    
    % Estimate image noise and mean brain intensity
    if strcmpi(modality,'MRI')        
        [sd(c),mu_brain] = spm_noise_estimate_mod(Nii_x(c),speak >= 2);
        
        mu(c)  = mu_brain;              % Mean brain intensity
        lam(c) = scl_lam/double(mu(c)); % This scaling is currently a bit arbitrary, and should be based on empiricism
    elseif strcmpi(modality,'CT')
        sd(c)  = noise_estimate_ct(Nii_x(c));
        
        mu(c)  = 0;
        lam(c) = lam_ct;
    end
    
    tau(c) = 1/(sd(c).^2); % Noise precision
end

% For decreasing regularisation with iteration number
lam0      = lam;
def       = spm_shoot_defaults;
sched_lam = def.sched;
sched_lam = sched_lam(end - min(numel(sched_lam) - 1,nit):end);

% lam = prod(vx_sr)*lam; % Scale regularisation with voxel size (only for super-resolution)
lam = sched_lam(1)*lam;

if rho == 0
    % Estimate rho (this value seems to lead to reasonably good convergence)
    rho = sqrt(mean(tau))/mean(lam);
end

if speak  >= 1
    % Print estimates
    fprintf('Estimated parameters are:\n');
    for c=1:C        
        fprintf('c=%i | sd=%f, mu=%f | tau=%f, lam=%f, rho=%f\n', c, sd(c), mu(c), tau(c), lam(c), rho);
    end
    fprintf('\n');
end

if coreg
    % Co-register input images
    Nii_x = coreg_ims(Nii_x,dir_tmp);
end

%--------------------------------------------------------------------------
% Allocate temporary variables
%--------------------------------------------------------------------------

if do_readwrite
    % Read/write temporary variables from disk (stored as NIfTIs)
    Nii_y = nifti;
    Nii_u = nifti;
    Nii_w = nifti;
else
    % Keep temporary variables in memory
    Nii_y = struct;
    Nii_u = struct;
    Nii_w = struct;
end
for c=1:C
    if do_readwrite
        fname_y = fullfile(dir_tmp,['y' num2str(c) '.nii']);
        fname_u = fullfile(dir_tmp,['u' num2str(c) '.nii']); 
        fname_w = fullfile(dir_tmp,['w' num2str(c) '.nii']);

        create_nii(fname_y,zeros(dm,    'single'),mat,[spm_type('float32') spm_platform('bigend')],'y');
        create_nii(fname_u,zeros([dm 3 2],'single'),mat,[spm_type('float32') spm_platform('bigend')],'u');
        create_nii(fname_w,zeros([dm 3 2],'single'),mat,[spm_type('float32') spm_platform('bigend')],'w');

        Nii_y(c) = nifti(fname_y);
        Nii_u(c) = nifti(fname_u);
        Nii_w(c) = nifti(fname_w);
    else
        Nii_y(c).dat = zeros(dm,    'single');
        Nii_u(c).dat = zeros([dm 3 2],'single');
        Nii_w(c).dat = zeros([dm 3 2],'single');
    end
end

%--------------------------------------------------------------------------
% Initialise variables
%--------------------------------------------------------------------------

infnrm = zeros(1,C);
dat    = zeros(1,C);
L      = [];
if strcmpi(method,'superres')   
    %---------------------------
    % Super-resolution
    %---------------------------
            
    % Initialise dat struct with projection matrices, etc.
    dat = init_dat(Nii_x,mat,dm);
            
    if superes_fmg
        % Compute infinity norm
        infnrm = zeros(1,C);
        for c=1:C        
            tmp       = At(A(ones(dat(c).dm,'single'),dat(c)),dat(c));
            infnrm(c) = max(tmp(:)); % If A is all positive, max(A'*A*ones(N,1)) gives the infinity norm
        end        
        clear tmp
    else
        % Define Laplace prior (in Fourier space)
        L = laplace_prior(dm,vx);
    end       
end

% Manage parfor
num_workers                        = min(C,num_workers);
if C == 1,             num_workers = 0; end
if num_workers == Inf, num_workers = nbr_parfor_workers; end
if num_workers > 1,    manage_parpool(num_workers);  end

msk = cell(1,C); % For saving locations of missing values so that they can be 're-applied' once the algorithm has finished
% for c=1:C
parfor (c=1:C,num_workers)  

    spm_field('boundary',1) % Set up boundary conditions that match the gradient operator
            
    % Observed image
    x = get_nii(Nii_x(c));  
        
    % Initial guesses for solution    
    if strcmpi(method,'superres')    
        %---------------------------
        % Super-resolution
        %---------------------------
        
        if superes_fmg
            y = get_nii(Nii_y(c));  
            for gnit=1:1 % Iterate Gauss-Newton

                % Gradient      
                Ayx = A(y,dat(c));
                for n=1:dat(c).N
                   Ayx{n} = Ayx{n} - x;
                end                
                rhs = At(Ayx,dat(c),tau(c))*(1/rho); 
                Ayx = [];
                rhs = rhs + spm_field('vel2mom',y,[vx 0 lam(c)^2 0]);

                % Hessian
                lhs = infnrm(c)*ones(dm,'single')*tau(c)/rho;

                % Compute GN step
                y   = y - spm_field(lhs,rhs,[vx 0 lam(c)^2 0 2 2]);
                lhs = [];
                rhs = [];

            end
            msk{c} = isfinite(y);
        else
            [y,msk{c}] = get_y_superres(Nii_x(c),dat(c),dm,mat); % 4th order b-splines
        end

    else  
        %---------------------------
        % Denoising
        %---------------------------        
    
        y = spm_field(tau(c)*ones(dm,'single'),tau(c)*x,[vx 0 lam(c)^2 0 2 2]); 
        
        if strcmpi(modality,'CT')
            msk{c} = isfinite(x) & x ~= 0 & x ~= min(x);
        else
            msk{c} = isfinite(x) & x ~= 0;
        end
    end        
    x = [];
    
    if strcmpi(modality,'MRI')
        % Ensure non-negativity
        y(y < 0) = 0;
    end 
        
    Nii_y(c) = put_nii(Nii_y(c),y);                
    y        = [];
    
    % Proximal variables
    Nii_u(c) = put_nii(Nii_u(c),zeros([dm 3 2],'single'));
    Nii_w(c) = put_nii(Nii_w(c),zeros([dm 3 2],'single'));
end
rhs = [];

%--------------------------------------------------------------------------
% Start solving
%--------------------------------------------------------------------------

if speak >= 1
    if tol == 0
        fprintf('Start %s, running %d iterations\n', method, nit);
    else
        fprintf('Start %s, running (max) %d iterations\n', method, nit);
    end
    tic; 
end

armijo = ones(1,C);
ll     = -Inf;
for it=1:nit
        
    % Decrease regularisation with iteration number
    lam = sched_lam(min(it,numel(sched_lam)))*lam0;    
    
    %------------------------------------------------------------------
    % Proximal operator for u
    % Here we solve for the MTV term of the objective function using
    % vectorial soft thresholding 
    %------------------------------------------------------------------

    unorm = 0;    
%     for c=1:C
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
    
    unorm = sqrt(unorm);
    scale = max(unorm - 1/rho,0)./(unorm + eps);
    clear unorm
        
    ll1 = zeros(1,C);
    ll2 = 0;
%     for c=1:C
    parfor (c=1:C,num_workers) % Loop over channels
    
        spm_field('boundary',1) % Set up boundary conditions that match the gradient operator
                
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
        %------------------------------------------------------------------
        rhs = [];
        if ~superes_fmg
            rhs = u - w/rho; 
            rhs = lam(c)*imdiv(rhs,vx);
        end
        
        x = get_nii(Nii_x(c)); % Get observed image
        
        if strcmpi(method,'superres')              
            %---------------------------
            % Super-resolution
            %---------------------------
                                  
            if superes_fmg
                %---------------------------
                % FMG
                % Here we want to solve for y using Gauss-Newton (GN)
                % optimisation
                %---------------------------      
                                        
                y = get_nii(Nii_y(c)); % Get solution
                         
                if armijo(c) < eps('single')
                    % At optimum
                    continue;
                end
                            
                for gnit=1:1 % Iterate Gauss-Newton
                    
                    % Gradient      
                    rhs       = w/rho - u; 
                    rhs       = lam(c)*imdiv(rhs,vx);
                    Ayx       = A(y,dat(c));
                    for n=1:dat(c).N
                       Ayx{n} = Ayx{n} - x;
                    end                
                    rhs       = rhs + At(Ayx,dat(c),tau(c))*(1/rho); 
                    Ayx       = [];
                    rhs       = rhs + spm_field('vel2mom',y,[vx 0 lam(c)^2 0]);

                    % Hessian
                    lhs = infnrm(c)*ones(dm,'single')*tau(c)/rho;

                    % Compute GN step
                    y   = y - spm_field(lhs,rhs,[vx 0 lam(c)^2 0 2 2]);
                    lhs = [];
                    rhs = [];
                end
            else
                %---------------------------
                % CG
                % Here we want to solve the linear system: lhs\rhs
                %---------------------------    
                
                % RHS
                rhs = rhs + At({x},dat(c),tau(c))*(1/rho); 
                     
                % LHS
                lhs = @(y) AtA(y,@(y) L(y),tau(c)/rho,lam(c)^2,dat(c));
                
                % Compute new y
                [y,it_cg,d_cg,t_cg] = cg_im_solver(lhs,rhs,get_nii(Nii_y(c)),nit_cg,tol_cg);

                if speak >= 1 
                    fprintf('%2d | %2d %10.1f %10.1f\n', c, it_cg, d_cg, t_cg); 
                end                        
            end
        else            
            %---------------------------
            % Denoising
            %---------------------------
            
            % RHS
            rhs = rhs + x*(tau(c)/rho);
            
            % LHS
            lhs = ones(dm,'single')*tau(c)/rho;
            
            % Compute new y
            y = spm_field(lhs,rhs,[vx 0 lam(c)^2 0 2 2]);
        end                               
        lhs = [];
        rhs = [];

        if strcmpi(modality,'MRI')
            % Ensure non-negativity
            y(y < 0) = 0;
        end 
        
        % Compute likelihood term of posterior
        ll1(c) = get_ll(method,y,x,tau(c),dat(c));
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
    clear scale
    
    % Compute log-posterior (objective value)
    ll2  = -sum(sum(sum(sqrt(ll2))));    
    ll   = [ll, sum(ll1) + ll2];    
    gain = abs((ll(end - 1)*(1 + 10*eps) - ll(end))/ll(end));
        
    % Some (potential) verbose               
    if speak >= 1, fprintf('%2d | %10.1f %10.1f %10.1f %0.6f\n', it, sum(ll1), ll2, sum(ll1) + ll2, gain); end
    if speak >= 2, show_progress(method,modality,ll,Nii_x,Nii_y,dm,nr,nc); end
    
    if tol > 0 && gain < tol && it > numel(sched_lam)
        % Finished
        break
    end
end

if speak >= 1, toc; end

%--------------------------------------------------------------------------
% Write results
%--------------------------------------------------------------------------

if strcmpi(method,'superres'), prefix = 'sr';
else,                          prefix = 'den';
end
   
Nii = nifti;
for c=1:C
    VO          = spm_vol(Nii_x(c).dat.fname);
    [~,nam,ext] = fileparts(VO.fname);
    VO.fname    = fullfile(dir_out,[prefix '_' nam ext]);
    VO.dim(1:3) = dm(1:3);    
    VO.mat      = mat;
    VO          = spm_create_vol(VO);
        
    Nii(c) = nifti(VO.fname);
    y      = get_nii(Nii_y(c));  
    if strcmpi(method,'superres')
        % Rescale intensities
        vx0 = sqrt(sum(Nii_x(c).mat(1:3,1:3).^2)); 
        scl = prod(vx0./vx);
        y   = scl*y;
    end
    if zeroMissing
        y(~msk{c}) = 0; % 'Re-apply' missing values        
    end
    Nii(c).dat.scl_slope = max(y(:))/1600;
    create(Nii(c));
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

if do_clean && (do_readwrite || (coreg && C > 1))
    % Clean-up temporary files
    rmdir(dir_tmp,'s');
end
%==========================================================================
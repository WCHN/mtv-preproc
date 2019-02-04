function Nii = spm_mtv_preproc(varargin)
% Multi-channel total variation (MTV) preprocessing of MR and CT data. 
%
% Requires that the SPM software is on the MATLAB path.
% SPM is available from: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
%
% For super-resolution, remember to compile private/pushpull.c (see 
% private/compile_pushpull)
%
% FORMAT Nii = spm_mtv_preproc(...)
%
% KEYWORD
% -------
%
% InputImages            - Cell array of either NIfTI filenames or nifti 
%                          objects. The cell array is of size 1 x C, where 
%                          C are the number of image channels. Each array 
%                          entry contains N_c images of the same channel. 
%                          If empty, uses spm_select ['']
% IterMax                - Maximum number of iteration 
%                          [method=superres:40, method=denoise:20]
% ADMMStepSize           - The infamous ADMM step size, set to zero for an 
%                          educated guess [0.1]
% Tolerance              - Convergence threshold, set to zero to run until 
%                          IterMax [0]
% RegScaleSuperResMRI    - Scaling of regularisation for MRI super-
%                          resolution [0.01]
% RegScaleDenoisingMRI    -Scaling of regularisation for MRI denoising [3.2]
% RegSuperresCT          - Regularisation used for CT super-resolution [0.001]
% RegDenoisingCT         - Regularisation used for CT denoising [0.06]
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
% VoxelSize              - Voxel size of super-resolved image. If empty, 
%                          sets voxel size to smallest voxel size in input 
%                          data [1]
% IterMaxCG              - Maximum number of iterations for conjugate gradient 
%                          solver used for super-resolution [12]
% ToleranceCG            - Convergence threshold for conjugate gradient 
%                          solver used for super-resolution [1e-3]
% CoRegister             - For super-resolution, co-register input images [true] 
% Modality               - Either MRI (denoise and super-resolution) or CT 
%                          (denoise) ['MRI']
% ReadWrite              - Keep variables in workspace (requires more RAM,
%                          but faster), or read/write from disk (requires 
%                          less RAM, but slower) [false] 
% ZeroMissingValues      - Set NaNs and zero values to zero after algorithm 
%                          has finished [C=1:true, C>1:false]
% IterGaussNewton        - Number of Gauss-Newton iterations for FMG 
%                          super-resolution [1]
% Reference              - Cell array (1xC) with reference images, if given
%                          computes PSNR and displays for each iteration of
%                          the algoirthm [{}]
% DecreasingReg          - Regularisation decreases over iterations, based
%                          on the scheduler in spm_shoot_defaults 
%                          [method=superres:true, method=denoise:false]
%
% OUTPUT
% ------
% 
% Nii - nifti object containing denoised/super-resolved images
% 
%__________________________________________________________________________
%
% Example: Super-resolve a set thick-sliced MRIs simulated from an IXI subject
%
% Simulate thick-sliced from IXI references by running the script:
% >> GenerateTestData % Down-sampling factor set by DownSampling parameter
%
% Read simulated thick-sliced IXI MRIs
% InputImages{1} = nifti(char({'./LowResData/ds_n1_IXI002-Guys-0828-PD.nii', ...
%                              './LowResData/ds_n2_IXI002-Guys-0828-PD.nii'}));
% InputImages{2} = nifti(char({'./LowResData/ds_n1_IXI002-Guys-0828-T2.nii', ...
%                              './LowResData/ds_n2_IXI002-Guys-0828-T2.nii'}));
% InputImages{3} = nifti(char({'./LowResData/ds_n1_IXI002-Guys-0828-T1.nii'}));
%
% Super-resolve the MRIs
% >> spm_mtv_preproc('InputImages',InputImages,'Method','superres','Verbose',2);
%
% Compare super-resolved with known ground-truth
% >> files_sr  = spm_select('FPList','./out', '^sr_.*\.nii$');
% >> files_ref = spm_select('FPList','./data','^IXI.*\.nii$');
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
spm_check_path;

%--------------------------------------------------------------------------
% Parse input
%--------------------------------------------------------------------------

p              = inputParser;
p.FunctionName = 'spm_mtv_preproc';
p.addParameter('InputImages', {}, @(in) ( isa(in,'nifti') || isempty(in) || ...
                                        ((ischar(in{1}) || isa(in{1},'nifti')) || (ischar(in{1}{1}) || isa(in{1}{1},'nifti'))) ) );
p.addParameter('IterMax', 0, @(in) (isnumeric(in) && in >= 0));
p.addParameter('ADMMStepSize', 0.1, @(in) (isnumeric(in) && in >= 0));
p.addParameter('Tolerance', 0, @(in) (isnumeric(in) && in >= 0));
p.addParameter('RegScaleSuperResMRI', 0.01, @(in) (isnumeric(in) && in > 0));
p.addParameter('RegScaleDenoisingMRI', 3.2, @(in) (isnumeric(in) && in > 0));
p.addParameter('RegSuperresCT', 0.001, @(in) (isnumeric(in) && in > 0));
p.addParameter('RegDenoisingCT', 0.06, @(in) (isnumeric(in) && in > 0));
p.addParameter('WorkersParfor', Inf, @(in) (isnumeric(in) && in >= 0));
p.addParameter('TemporaryDirectory', 'tmp', @ischar);
p.addParameter('OutputDirectory', 'out', @ischar);
p.addParameter('Method', 'denoise', @(in) (ischar(in) && (strcmpi(in,'denoise') || strcmpi(in,'superres'))));
p.addParameter('Verbose', 1, @(in) (isnumeric(in) && in >= 0 && in <= 3));
p.addParameter('CleanUp', true, @islogical);
p.addParameter('VoxelSize', [], @(in) ((isnumeric(in) && (numel(in) == 1 || numel(in) == 3)) && ~any(in <= 0)) || isempty(in));
p.addParameter('IterMaxCG', 12, @(in) (isnumeric(in) && in > 0));
p.addParameter('ToleranceCG', 1e-4, @(in) (isnumeric(in) && in >= 0));
p.addParameter('CoRegister', true, @islogical);
p.addParameter('Modality', 'MRI', @(in) (ischar(in) && (strcmpi(in,'MRI') || strcmpi(in,'CT'))));
p.addParameter('ReadWrite', false, @islogical);
p.addParameter('ZeroMissingValues', [], @(in) (islogical(in) || isnumeric(in)));
p.addParameter('IterGaussNewton', 1, @(in) (isnumeric(in) && in > 0));
p.addParameter('Reference', {}, @iscell);
p.addParameter('DecreasingReg', [], @(in) (islogical(in) || isempty(in)));
p.parse(varargin{:});
InputImages  = p.Results.InputImages;
nit          = p.Results.IterMax;
tol          = p.Results.Tolerance;
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
do_readwrite = p.Results.ReadWrite; 
rho          = p.Results.ADMMStepSize; 
zeroMissing  = p.Results.ZeroMissingValues; 
nitgn        = p.Results.IterGaussNewton; 
ref          = p.Results.Reference; 
dec_reg      = p.Results.DecreasingReg; 

%--------------------------------------------------------------------------
% Preliminaries
%--------------------------------------------------------------------------

if ~isempty(ref) && strcmpi(method,'superres')
    error('Super-resolution with reference image(s) not yet implemented!');
end

% Get image data
[Nii_x,C,N0] = parse_input_data(InputImages,method);

if isempty(zeroMissing)    
    % Missing values (NaNs and zeros) will be...
    if C == 1 && numel(Nii_x{1}) == 1
        % ...set to zero after algorithm finishes
        zeroMissing = true;
    else
        % ...filled in by the algorithm
        zeroMissing = false;
    end
end

% Super-resolution voxel-size related
if isempty(vx_sr)
    % Get voxel-size from input images
    vx_sr = get_vx_sr(Nii_x);
elseif numel(vx_sr) == 1
    vx_sr = vx_sr*ones(1,3); 
end
if vx_sr(1) < 0.9
    % Voxels are quite small, read-write aux. variables to not run in to
    % memory issues..
    do_readwrite = true;
end

% Make some directories
if  exist(dir_tmp,'dir') == 7,  rmdir(dir_tmp,'s'); end
if  do_readwrite || (coreg && (C > 1 || numel(Nii_x{1}) > 1)), mkdir(dir_tmp); end
if ~(exist(dir_out,'dir') == 7),  mkdir(dir_out);  end

% Set defaults, get voxel size, orientation matrix and image dimensions
if strcmpi(method,'denoise')
    %---------------------------
    % Denoising
    %---------------------------
    
    if isempty(dec_reg), dec_reg = false; end
    if nit == 0,         nit     = 20; end
    
    scl_lam = p.Results.RegScaleDenoisingMRI;
    lam_ct  = p.Results.RegDenoisingCT; 
    
    mat = Nii_x{1}(1).mat;
    vx  = sqrt(sum(mat(1:3,1:3).^2));   
    dm  = Nii_x{1}(1).dat.dim;
elseif strcmpi(method,'superres')
    %---------------------------
    % Super-resolution
    %---------------------------
            
    if isempty(dec_reg), dec_reg = true; end
    if nit == 0,         nit     = 40; end
    
    scl_lam = p.Results.RegScaleSuperResMRI;
    lam_ct  = p.Results.RegSuperresCT;
    
    % For super-resolution, calculate orientation matrix and dimensions 
    % from maximum bounding-box
    vx       = vx_sr;
    [mat,dm] = max_bb_orient(Nii_x,vx);
end

%--------------------------------------------------------------------------
% Estimate model hyper-parameters
%--------------------------------------------------------------------------

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
nr = floor(sqrt(N0));
nc = ceil(N0/nr);  
    
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
def       = spm_shoot_defaults;
sched_lam = def.sched;
sched_lam = sched_lam(1:min(numel(sched_lam),nit));
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
    clear atau
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

%--------------------------------------------------------------------------
% Co-register input images
%--------------------------------------------------------------------------

if coreg
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
            
    % Compute infinity norm
    infnrm = zeros(1,C);
    for c=1:C        
        tmp       = At(A(ones(dat(c).dm,'single'),dat(c)),dat(c));
        infnrm(c) = max(tmp(:)); % If A is all positive, max(A'*A*ones(N,1)) gives the infinity norm
    end        
    clear tmp 
end

% Manage parfor
num_workers                        = min(C,num_workers);
if C == 1,             num_workers = 0; end
if num_workers == Inf, num_workers = nbr_parfor_workers; end
if num_workers > 1,    manage_parpool(num_workers);  end

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

ll = -Inf;
for it=1:nit
        
    if dec_reg
        % Decrease regularisation with iteration number
        lam = sched_lam(min(it,numel(sched_lam)))*lam0;    
    end
    
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
    
    unorm = sqrt(unorm);
    scale = max(unorm - 1/rho,0)./(unorm + eps);
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
        
        u = bsxfun(@times,u,scale);
            
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
                    mskn          = x{n} ~= 0 & isfinite(x{n});
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
    clear scale
    
    % Compute log-posterior (objective value)
    ll2  = -sum(sum(sum(sqrt(ll2))));    
    ll   = [ll, sum(ll1) + ll2];    
    gain = abs((ll(end - 1)*(1 + 10*eps) - ll(end))/ll(end));
        
    % Some (potential) verbose               
    if speak >= 1 || ~isempty(ref)
        if ~isempty(ref)
            % Reference image(s) given, compute and output PSNR
            psnrs = zeros(1,C);
            for c=1:C
                psnrs(c) = get_psnr(get_nii(Nii_y(c)),ref{c});
            end
            
            fprintf('%2d | %10.1f %10.1f %10.1f %0.6f|%s\n', it, sum(ll1), ll2, sum(ll1) + ll2, gain, sprintf(' %2.2f', psnrs)); 
        else
            fprintf('%2d | %10.1f %10.1f %10.1f %0.6f\n', it, sum(ll1), ll2, sum(ll1) + ll2, gain); 
        end
        
        if speak >= 2
            show_progress(method,modality,ll,Nii_x,Nii_y,dm); 
        end
    end   
    
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
    % Set output filename
    [~,nam,ext] = fileparts(Nii_x{c}(1).dat.fname);
    nfname      = fullfile(dir_out,[prefix '_' nam ext]);
    
    % Get output image data
    y = get_nii(Nii_y(c));  

    if zeroMissing
        y(~msk{c}) = 0; % 'Re-apply' missing values        
    end
    
    % Write to NIfTI
    Nii(c) = create_nii(nfname,y,mat,[spm_type('float32') spm_platform('bigend')],'MTV recovered');
end

%--------------------------------------------------------------------------
% Show input and solved
%--------------------------------------------------------------------------

if speak >= 3
    fnames = cell(1,2*C);
    cnt    = 1;
    for c=1:2:2*C    
        fnames{c}     = Nii_x0{cnt}(1).dat.fname;    
        fnames{c + 1} = Nii(cnt).dat.fname;
        cnt           = cnt + 1;
    end

    spm_check_registration(char(fnames))
end

if do_clean && (do_readwrite || (coreg && (C > 1 || numel(Nii_x{1}) > 1)))
    % Clean-up temporary files
    rmdir(dir_tmp,'s');
end
%==========================================================================
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
% InputImages          - Cell array of either NIfTI filenames or nifti 
%                        objects. The cell array is of size 1 x C, where 
%                        C are the number of image channels. Each array 
%                        entry contains N_c images of the same channel. 
%                        If empty, uses spm_select ['']
% IterMax              - Maximum number of iteration 
%                        [method=superres:40, method=denoise:20]
% ADMMStepSize         - The infamous ADMM step size, set to zero for an 
%                        educated guess [0]
% Tolerance            - Convergence threshold, set to zero to run until 
%                        IterMax [0]
% RegScaleSuperResMRI  - Scaling of regularisation for MRI super-
%                        resolution [20]
% RegScaleDenoisingMRI - Scaling of regularisation for MRI denoising, 
%                        increase this value for stronger denoising [3.2]
% WorkersParfor        - Maximum number of parfor workers [Inf]
% TemporaryDirectory   - Directory for temporary files ['./tmp']
% OutputDirectory      - Directory for denoised images ['./out']
% Method               - Does either denoising ('denoise') or 
%                        super-resolution ('superres') ['denoise']
% Verbose              - Verbosity level: 
%                        *  0  = quiet
%                        * [1] = write   (print objective value and parameter estimates)
%                        *  2  = draw   +(figure w. log likelihood, mixture fits, recons)
%                        *  3  = result +(show observed and reconstructed images 
%                                         in spm_check_registration, when finished)
% CleanUp              - Delete temporary files [true] 
% VoxelSize            - Voxel size of super-resolved image [1 1 1]
% CoRegister           - For super-resolution, co-register input images [true] 
% Modality             - Either MRI (denoise and super-resolution) or CT 
%                        (denoise) ['MRI']
% RegSuperresCT        - Regularisation used for CT denoising [0.001]
% RegDenoisingCT       - Regularisation used for CT super-resolution [0.04]
% ReadWrite            - Keep variables in workspace (requires more RAM,
%                        but faster), or read/write from disk (requires 
%                        less RAM, but slower) [false] 
% ZeroMissingValues    - Set NaNs and zero values to zero after algorithm 
%                        has finished [C=1:true, C>1:false]
% IterGaussNewton      - Number of Gauss-Newton iterations for FMG 
%                        super-resolution [1]
% Reference            - Struct with NIfTI reference images, if given 
%                        computes PSNR and displays it, for each iteration of 
%                        the algoirthm [{}] 
% DecreasingReg        - Regularisation decreases over iterations, based
%                        on the scheduler in spm_shoot_defaults 
%                        [method=superres:true, method=denoise:false]
% SliceProfile         - Slice selection profile along directions x, y, z.
%                        It should have the same shape as InputImages.
%                        *  1  = Gaussian  (FWHM = low-resolution voxel)
%                        * [2] = Rectangle (length = low-resolution voxel)
% SliceGap             - Gap between slices in mm along directions x, y, z. 
%                        A positive value means a gap, a negative value 
%                        means an overlap. [0]
% EstimateRigid        - Optimise a rigid alignment between observed images
%                        and their corresponding channel's reconstruction
%                        [false]
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
spm_check_path('Longitudinal');

%--------------------------------------------------------------------------
% Parse input
%--------------------------------------------------------------------------

p              = inputParser;
p.FunctionName = 'spm_mtv_preproc';
p.addParameter('InputImages', {}, @(in) ( isa(in,'nifti') || isempty(in) || ...
                                        ((ischar(in{1}) || isa(in{1},'nifti')) || (ischar(in{1}{1}) || isa(in{1}{1},'nifti'))) ) );
p.addParameter('IterMax', 0, @(in) (isnumeric(in) && in >= 0));
p.addParameter('ADMMStepSize', 0, @(in) (isnumeric(in) && in >= 0));
p.addParameter('Tolerance', 0, @(in) (isnumeric(in) && in >= 0));
p.addParameter('RegScaleSuperResMRI', 20, @(in) (isnumeric(in) && in > 0));
p.addParameter('RegScaleDenoisingMRI', 3.2, @(in) (isnumeric(in) && in > 0));
p.addParameter('RegSuperresCT', 0.001, @(in) (isnumeric(in) && in > 0));
p.addParameter('RegDenoisingCT', 0.06, @(in) (isnumeric(in) && in > 0));
p.addParameter('WorkersParfor', Inf, @(in) (isnumeric(in) && in >= 0));
p.addParameter('TemporaryDirectory', 'tmp', @ischar);
p.addParameter('OutputDirectory', 'out', @ischar);
p.addParameter('Method', 'denoise', @(in) (ischar(in) && (strcmpi(in,'denoise') || strcmpi(in,'superres'))));
p.addParameter('Verbose', 1, @(in) (isnumeric(in) && in >= 0 && in <= 3));
p.addParameter('CleanUp', true, @islogical);
p.addParameter('VoxelSize', [1 1 1], @(in) ((isnumeric(in) && (numel(in) == 1 || numel(in) == 3)) && ~any(in <= 0)) || isempty(in));
p.addParameter('CoRegister', true, @islogical);
p.addParameter('Modality', 'MRI', @(in) (ischar(in) && (strcmpi(in,'MRI') || strcmpi(in,'CT'))));
p.addParameter('ReadWrite', false, @islogical);
p.addParameter('ZeroMissingValues', [], @(in) (islogical(in) || isnumeric(in)));
p.addParameter('IterGaussNewton', 1, @(in) (isnumeric(in) && in > 0));
p.addParameter('Reference', [], @(in)  isa(in,'nifti'));
p.addParameter('DecreasingReg', [], @(in) (islogical(in) || isempty(in)));
p.addParameter('SliceProfile', 2, @(in) (isnumeric(in) || iscell(in)));
p.addParameter('SliceGap', 0, @(in) (isnumeric(in) || iscell(in)));
p.addParameter('EstimateRigid', false, @islogical);
p.parse(varargin{:});
InputImages   = p.Results.InputImages;
nit           = p.Results.IterMax;
tol           = p.Results.Tolerance;
num_workers   = p.Results.WorkersParfor;
dir_tmp       = p.Results.TemporaryDirectory;
dir_out       = p.Results.OutputDirectory;
method        = p.Results.Method;
speak         = p.Results.Verbose; 
do_clean      = p.Results.CleanUp; 
vx_sr         = p.Results.VoxelSize; 
coreg         = p.Results.CoRegister; 
modality      = p.Results.Modality; 
do_readwrite  = p.Results.ReadWrite; 
zeroMissing   = p.Results.ZeroMissingValues; 
Nii_ref       = p.Results.Reference; 
dec_reg       = p.Results.DecreasingReg;
window        = p.Results.SliceProfile;
gap           = p.Results.SliceGap;
EstimateRigid = p.Results.EstimateRigid;

%--------------------------------------------------------------------------
% Preliminaries
%--------------------------------------------------------------------------

% Some sanity checks
if ~isempty(Nii_ref) && strcmpi(method,'superres')
    error('Super-resolution with reference image(s) not yet implemented!');
end
if EstimateRigid && strcmpi(method,'denoise')
    error('Denoising with rigid alignment not yet implemented!');
end

% Get image data
[Nii_x,C] = parse_input_data(InputImages,method);

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

% Manage parfor
num_workers                        = min(C,num_workers);
if C == 1,             num_workers = 0; end
if num_workers == Inf, num_workers = nbr_parfor_workers; end
if num_workers > 1,    manage_parpool(num_workers);  end
    
%--------------------------------------------------------------------------
% Co-register input images (modifies images' orientation matrices)
%--------------------------------------------------------------------------

if coreg
    Nii_x = coreg_ims(Nii_x,dir_tmp);
end

%--------------------------------------------------------------------------
% Initialise super-resolution/denoising
%--------------------------------------------------------------------------

% Set defaults, such as, voxel size, orientation matrix and image dimensions
if strcmpi(method,'denoise')
    %---------------------------
    % Denoising
    %---------------------------
    
    if isempty(dec_reg), dec_reg = false; end
    if nit == 0,         nit     = 20; end
    
    mat    = Nii_x{1}(1).mat;
    vx     = sqrt(sum(mat(1:3,1:3).^2));   
    dm     = Nii_x{1}(1).dat.dim;
    dat    = struct;
    infnrm = 0;
elseif strcmpi(method,'superres')
    %---------------------------
    % Super-resolution
    %---------------------------
            
    if isempty(dec_reg), dec_reg = true; end
    if nit == 0,         nit     = 40; end
    
    % For super-resolution, calculate orientation matrix and dimensions 
    % from maximum bounding-box
    vx       = vx_sr;
    [mat,dm] = max_bb_orient(Nii_x,vx);
    
    % Initialise dat struct with projection matrices, etc.
    dat = init_dat(Nii_x,mat,dm,window,gap);            
end

%--------------------------------------------------------------------------
% Estimate model hyper-parameters
%--------------------------------------------------------------------------

[tau,lam,rho,sched_lam,lam0,Nii_x0] = estimate_model_hyperpars(Nii_x,dec_reg,nit,p);

%--------------------------------------------------------------------------
% Allocate temporary variables
%--------------------------------------------------------------------------

[Nii_y,Nii_u,Nii_w,Nii_H] = alloc_aux_vars(do_readwrite,C,dm,mat,dir_tmp);

if strcmpi(method,'superres')
    % Compute approximation to the diagonal of the Hessian 
    for c=1:C        
        H        = At(A(ones(dat(c).dm,'single'),dat(c)),dat(c));        
        Nii_H(c) = put_nii(Nii_H(c),H);
    end            
    clear H
end

%--------------------------------------------------------------------------
% Create intial estimate of solution (y)
%--------------------------------------------------------------------------

[Nii_y,msk] = estimate_initial_y(Nii_x,Nii_y,Nii_H,dat,tau,rho,lam,vx,dm,num_workers,p);

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

% For storing rigid optimisation line-search parameters
armijo_rigid = cell(1,C);
for c=1:C
    N               = numel(Nii_x{c});
    armijo_rigid{c} = ones([1 N]);
end

% For storing model log-likelihood, for each iteration
ll = -Inf; 

for it=1:nit % Start iterating
        
    if dec_reg
        % Decrease regularisation with iteration number
        lam = sched_lam(min(it,numel(sched_lam)))*lam0;    
    end
    
    %----------------------------------------------------------------------
    % Update recovered image(s) (Nii_y)
    %----------------------------------------------------------------------
    
    [Nii_y,Nii_u,Nii_w,ll1,ll2]= update_y(Nii_x,Nii_y,Nii_u,Nii_w,Nii_H,dat,tau,rho,lam,vx,dm,num_workers,p);
   
    % Compute log-posterior (objective value)        
    ll   = [ll, -(sum(ll1) + ll2)]; % Minus sign because YB wants to see increasing objective functions...
    gain = get_gain(ll);
                       
    if speak >= 1 || ~isempty(Nii_ref)
        % Some verbose    
        
        if ~isempty(Nii_ref)
            % Reference image(s) given, compute and output PSNR
            psnrs = zeros(1,C);
            ssims = zeros(1,C);
            for c=1:C
                psnrs(c) = get_psnr(get_nii(Nii_y(c)),get_nii(Nii_ref(c)));
                ssims(c) = ssim(get_nii(Nii_y(c)),get_nii(Nii_ref(c)));
            end
            
            fprintf('%2d | ll1=%10.1f, ll2=%10.1f, ll=%10.1f, gain=%0.6f | psnr =%s, ssim =%s\n', it, sum(ll1), ll2, sum(ll1) + ll2, gain, sprintf(' %2.2f', psnrs), sprintf(' %2.2f', ssims)); 
        else
            fprintf('%2d | ll1=%10.1f, ll2=%10.1f, ll=%10.1f, gain=%0.6f\n', it, sum(ll1), ll2, sum(ll1) + ll2, gain); 
        end
        
        if speak >= 2
            show_progress(method,modality,ll,Nii_x,Nii_y,dm); 
        end
    end   
    
    % Check convergence
    if tol > 0 && gain < tol && it > numel(sched_lam)
        % Finished
        break
    end
    
    if EstimateRigid        
        
        %------------------------------------------------------------------
        % Update rigid alignment matrices (dat(c).A(n).q)
        %------------------------------------------------------------------
                
        [dat,ll1,armijo_rigid] = update_rigid(Nii_x,Nii_y,dat,tau,armijo_rigid,num_workers,speak);
        
        % Compute log-posterior (objective value)        
        ll   = [ll, -(sum(ll1) + ll2)]; % Minus sign because YB wants to see increasing objective functions...    
        gain = get_gain(ll);
    
        if speak >= 1
            % Some verbose    
            fprintf('   | ll1=%10.1f, ll2=%10.1f, ll=%10.1f, gain=%0.6f\n', sum(ll1), ll2, sum(ll1) + ll2, gain); 
            
            if speak >= 2
                show_progress(method,modality,ll,Nii_x,Nii_y,dm); 
            end
        end                           
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
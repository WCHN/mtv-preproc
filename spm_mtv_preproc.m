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
% IterMax              - Maximum number of iteration [30]
% IterImage            - Maximum number of iterations for solving for the
%                        super-resolved image(s) [3]
% ADMMStepSize         - The infamous ADMM step size, set to zero for an 
%                        educated guess [0]
% Tolerance            - Convergence threshold, set to zero to run until 
%                        IterMax [1e-4]
% RegScaleSuperResMRI  - Scaling of regularisation for MRI super-
%                        resolution [6]
% RegScaleDenoisingMRI - Scaling of regularisation for MRI denoising, 
%                        increase this value for stronger denoising [6]
% RegSuperresCT        - Regularisation used for CT denoising [0.06]
% RegDenoisingCT       - Regularisation used for CT super-resolution [0.06]
% WorkersParfor        - Maximum number of parfor workers [Inf]
% TemporaryDirectory   - Directory for temporary files ['./tmp']
% OutputDirectory      - Directory for denoised images. If none given,
%                        writes results to same folder as input (prefixed 
%                        either den or sr) ['']
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
% ReadWrite            - Keep variables in workspace (requires more RAM,
%                        but faster), or read/write from disk (requires 
%                        less RAM, but slower) [false] 
% ZeroMissingValues    - Set NaNs and zero values to zero after algorithm 
%                        has finished [C=1:true, C>1:false]
% IterGaussNewtonImage - Number of Gauss-Newton iterations for solving for 
%                        super-resolution [1]
% IterGaussNewtonRigid - Number of Gauss-Newton iterations for solving for
%                        rigid registration [1]
% Reference            - Struct with NIfTI reference images, if given 
%                        computes PSNR and displays it, for each iteration of 
%                        the algoirthm [{}] 
% DecreasingReg        - Regularisation decreases over iterations, based
%                        on the scheduler in spm_shoot_defaults 
%                        [method=superres:true, method=denoise:false]
% SliceProfile         - Slice selection profile, either a scalar (and the 
%                        slice-select direction will be found automatically) 
%                        or along directions x, y, z (same shape as
%                        InputImages):
%                        * 1 = Gaussian  (FWHM   = low-resolution voxel)
%                        * 2 = Rectangle (length = low-resolution voxel)
%                        [In-plane: Gaussian, Through-plane: Rectangle]
% SliceGap             - Gap between slices, either a scalar (and the 
%                        slice-select direction will be found automatically) 
%                        or along directions x, y, z (same shape as
%                        InputImages). A positive value means a gap, a 
%                        negative value means an overlap. [0]
% SliceGapUnit         - Percentage ('%') or milimeters ('mm') ['%']
% EstimateRigid        - Optimise a rigid alignment between observed images
%                        and their corresponding channel's reconstruction
%                        [false]
% MeanCorrectRigid     - Mean correct the rigid-body transform parameters 
%                        q [false]
% PaddingBB            - Pad bounding box with extra zeros in each
%                        direction [0]
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

set_boundary_conditions;
    
%--------------------------------------------------------------------------
% Parse input
%--------------------------------------------------------------------------

p              = inputParser;
p.FunctionName = 'spm_mtv_preproc';
p.addParameter('InputImages', {}, @(in) ( isa(in,'nifti') || isempty(in) || ...
                                        ((ischar(in{1}) || isa(in{1},'nifti')) || (ischar(in{1}{1}) || isa(in{1}{1},'nifti'))) ) );
p.addParameter('IterMax', 30, @(in) (isnumeric(in) && in >= 0));
p.addParameter('IterImage', 3, @(in) (isnumeric(in) && in > 0));
p.addParameter('ADMMStepSize', 0, @(in) (isnumeric(in) && in >= 0));
p.addParameter('Tolerance', 1e-4, @(in) (isnumeric(in) && in >= 0));
p.addParameter('RegScaleSuperResMRI', 6, @(in) (isnumeric(in) && in > 0));
p.addParameter('RegScaleDenoisingMRI', 6, @(in) (isnumeric(in) && in > 0));
p.addParameter('RegSuperresCT', 0.06, @(in) (isnumeric(in) && in > 0));
p.addParameter('RegDenoisingCT', 0.06, @(in) (isnumeric(in) && in > 0));
p.addParameter('WorkersParfor', Inf, @(in) (isnumeric(in) && in >= 0));
p.addParameter('TemporaryDirectory', 'Temp', @ischar);
p.addParameter('OutputDirectory', '', @ischar);
p.addParameter('Method', 'denoise', @(in) (ischar(in) && (strcmpi(in,'denoise') || strcmpi(in,'superres'))));
p.addParameter('Verbose', 1, @(in) (isnumeric(in) && in >= 0 && in <= 3));
p.addParameter('CleanUp', true, @islogical);
p.addParameter('VoxelSize', [1 1 1], @(in) ((isnumeric(in) && (numel(in) == 1 || numel(in) == 3)) && ~any(in <= 0)) || isempty(in));
p.addParameter('CoRegister', true, @islogical);
p.addParameter('Modality', 'MRI', @(in) (ischar(in) && (strcmpi(in,'MRI') || strcmpi(in,'CT'))));
p.addParameter('ReadWrite', false, @islogical);
p.addParameter('ZeroMissingValues', [], @(in) (islogical(in) || isnumeric(in)));
p.addParameter('IterGaussNewtonImage', 1, @(in) (isnumeric(in) && in > 0));
p.addParameter('IterGaussNewtonRigid', 1, @(in) (isnumeric(in) && in > 0));
p.addParameter('Reference', {}, @(in)  (isa(in,'nifti') || isempty(in)));
p.addParameter('DecreasingReg', [], @(in) (islogical(in) || isempty(in)));
p.addParameter('SliceProfile', {}, @(in) (isnumeric(in) || iscell(in)));
p.addParameter('SliceGap', 0, @(in) (isnumeric(in) || iscell(in)));
p.addParameter('SliceGapUnit', '%', @(in) (ischar(in) && (strcmp(in,'%') || strcmp(in,'mm'))));
p.addParameter('EstimateRigid', false, @islogical);
p.addParameter('MeanCorrectRigid', true, @islogical);
p.addParameter('PaddingBB', 0, @isnumeric);
p.parse(varargin{:});
InputImages   = p.Results.InputImages;
nit           = p.Results.IterMax;
nity          = p.Results.IterImage;
tol           = p.Results.Tolerance;
num_workers   = p.Results.WorkersParfor;
dir_tmp       = p.Results.TemporaryDirectory;
dir_out       = p.Results.OutputDirectory;
method        = p.Results.Method;
speak         = p.Results.Verbose; 
do_clean      = p.Results.CleanUp; 
vx_sr         = p.Results.VoxelSize; 
coreg         = p.Results.CoRegister; 
do_readwrite  = p.Results.ReadWrite; 
zeroMissing   = p.Results.ZeroMissingValues; 
Nii_ref       = p.Results.Reference; 
dec_reg       = p.Results.DecreasingReg;
window        = p.Results.SliceProfile;
gap           = p.Results.SliceGap;
gapunit       = p.Results.SliceGapUnit;
EstimateRigid = p.Results.EstimateRigid;
bb_padding    = p.Results.PaddingBB;

%--------------------------------------------------------------------------
% Preliminaries
%--------------------------------------------------------------------------

% Struct that will hold model variables
Nii = struct;

% Get image data
[Nii,C,is3d] = parse_input_data(Nii,InputImages,method);

if isempty(zeroMissing)    
    % Missing values (NaNs and zeros) will be...
    if C == 1 && numel(Nii.x{1}) == 1
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
    vx_sr = get_vx_sr(Nii.x);
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
if  do_readwrite || (coreg && (C > 1 || numel(Nii.x{1}) > 1)), mkdir(dir_tmp); end
if ~(exist(dir_out,'dir') == 7),  mkdir(dir_out);  end

% Manage parfor
num_workers                        = min(C,num_workers);
if C == 1,             num_workers = 0; end
if num_workers == Inf, num_workers = nbr_parfor_workers; end
if num_workers > 1,    manage_parpool(num_workers);  end

if speak >= 3
    % So that Verbose = 3 works for superres (because Nii_x are copied, then copies are deleted)
    Nii.x0 = Nii.x; 
end
    
if coreg
    % Make copies input data and update Nii_x        
    Nii.x = copy_ims(Nii.x,dir_tmp);
end

% Flag saying if we solve using projection matrices (A, At), or not
use_projmat = ~(strcmpi(method,'denoise') && ~EstimateRigid);

% Some sanity checks
if ~isempty(Nii_ref) && use_projmat
    error('Solving with projection matrices and reference image(s) not yet implemented!');
end

%--------------------------------------------------------------------------
% Co-register input images (modifies images' orientation matrices)
%--------------------------------------------------------------------------

if coreg && is3d
    Nii.x = coreg_ims(Nii.x);
end

%--------------------------------------------------------------------------
% Initialise super-resolution/denoising
%--------------------------------------------------------------------------

% Set defaults, such as, voxel size, orientation matrix and image dimensions
if strcmpi(method,'denoise')
    %---------------------------
    % Denoising
    %---------------------------
    
    if isempty(dec_reg) && EstimateRigid, dec_reg = true; end
    if isempty(dec_reg),                  dec_reg = false; end
    
    vx = sqrt(sum(Nii.x{1}(1).mat(1:3,1:3).^2));
elseif strcmpi(method,'superres')
    %---------------------------
    % Super-resolution
    %---------------------------
            
    if isempty(dec_reg), dec_reg = true; end
    
    % For super-resolution, calculate orientation matrix and dimensions 
    % from maximum bounding-box
    vx = vx_sr;    
end

% Get recovered images' dimensions and otientation matrices
if use_projmat
    [mat,dm]  = max_bb_orient(Nii.x,vx,bb_padding);    
else
    mat       = Nii.x{1}(1).mat;
    dm        = Nii.x{1}(1).dat.dim;
    if ~is3d
        dm(3) = 1;
    end
end

% Initialise dat struct with projection matrices, etc.
dat = init_dat(method,Nii.x,mat,dm,window,gap,gapunit);    

% Allocate auxiliary variables
Nii = alloc_aux_vars(Nii,do_readwrite,C,dm,mat,dir_tmp);

if use_projmat
    % Compute approximation to the diagonal of the Hessian 
    Nii.H = approx_hessian(Nii.H,dat);
end

% Rigid optimisation line-search parameter
armijo = get_armijo(dat);

%--------------------------------------------------------------------------
% Estimate model hyper-parameters
%--------------------------------------------------------------------------

[tau,lam,rho,sched_lam,lam0] = estimate_model_hyperpars(Nii.x,dec_reg,nit,p);

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

if ~isempty(Nii_ref)
    % Reference image(s) given, compute SSIM and PSNR

    % Observed and reference
    [psnr1,ssim1] = compute_image_metrics(Nii.x,Nii_ref);
    fprintf('   | ll1=%10.1f, ll2=%10.1f, ll=%10.1f, gain=%0.6f | psnr=%2.3f, ssim=%1.3f\n', 0, 0, 0, 0, psnr1, ssim1); 

    % Initial solution and reference
    [psnr1,ssim1] = compute_image_metrics(Nii.y,Nii_ref);
    fprintf('%2d | ll=%10.1f, ll1=%10.1f, ll2=%10.1f, gain=%0.6f | psnr=%2.3f, ssim=%1.3f\n', 0, 0, 0, 0, psnr1, ssim1); 
end

ll = -Inf;
llpart = 1;
for it=1:nit % Start main loop
        
    if dec_reg
        % Decrease regularisation with iteration number
        if it > 1
            osched = sched_lam(min(it - 1,numel(sched_lam)));
        end
        sched  = sched_lam(min(it,numel(sched_lam)));
        lam    = sched*lam0;    
        if it > 1 && osched ~= sched
            % Regularisation changed, reset rigid armijo to one
            armijo = get_armijo(dat);             
        end
    else
        sched = 1;
    end
    
    %----------------------------------------------------------------------
    % ADMM to update image
    %----------------------------------------------------------------------    
    
    for ity=1:nity % Start y loop
        
        % Update Nii_y, Nii_w, Nii_u
        [Nii,ll1,ll2]= update_image(Nii,dat,tau,rho,lam,num_workers,p);

        % Compute log-posterior (objective value)        
        ll   = [ll, sum(ll1) + ll2];
        llpart = [llpart 1];
        gain = get_gain(ll);

        if speak >= 1 || ~isempty(Nii_ref)
            % Some verbose    

            if ~isempty(Nii_ref)
                % Reference image(s) given, compute SSIM and PSNR
                [psnr1,ssim1] = compute_image_metrics(Nii.y,Nii_ref);
                
                fprintf('%2d | ll=%10.1f, ll1=%10.1f, ll2=%10.1f, gain=%0.6f | psnr=%2.3f, ssim=%1.3f\n', it, ll(end), sum(ll1), ll2, gain, psnr1, ssim1); 
            else
                fprintf('%2d | ll=%10.1f, ll1=%10.1f, ll2=%10.1f, gain=%0.6f\n', it, ll(end), sum(ll1), ll2, gain); 
            end

            if speak >= 2
                show_model('ll',ll,llpart);                
            end
        end   
        
        if  tol > 0 && ity > 1 && ...
           (gain < 1e-3 && sched ~= 1 || ...
            gain < tol  && sched == 1)
            break
        end
        
    end % End y loop
    
    % Check convergence
    if tol > 0 && gain < tol && it > 1
        % Finished!
        break
    end
    
    if EstimateRigid
        
        %------------------------------------------------------------------
        % Gauss-Newton to update rigid alignment
        %------------------------------------------------------------------
                
        % Update q
        [dat,armijo,ll1] = update_rigid(Nii,dat,tau,armijo,num_workers,p);                            
        ll = [ll, sum(ll1) + ll2];
        llpart = [llpart 2];
        if speak >= 2
            show_model('ll',ll,llpart);                
        end
        
        % Update approximation to the diagonal of the Hessian 
        Nii.H        = approx_hessian(Nii.H,dat);
    end
 
end % End main loop

if speak >= 1, toc; end

%--------------------------------------------------------------------------
% Write results
%--------------------------------------------------------------------------

if strcmpi(method,'superres'), prefix = 'sr';
else,                          prefix = 'den';
end
   
Nii_out = nifti;
for c=1:C
    % Set output filename
    [pth,nam,ext] = fileparts(Nii.x{c}(1).dat.fname);
    if isempty(dir_out)
        nfname    = fullfile(pth,[prefix '_' nam ext]);
    else
        nfname    = fullfile(dir_out,[prefix '_' nam ext]);
    end
    
    % Get output image data
    y = get_nii(Nii.y(c));  
    
    % Write to NIfTI
    if ~is3d
        mat(3,4) = Nii.x{c}(1).mat(3,4);
    end
    Nii_out(c) = create_nii(nfname,y,mat,[spm_type('float32') spm_platform('bigend')],'MTV recovered');
end

%--------------------------------------------------------------------------
% Show input and solved
%--------------------------------------------------------------------------

if speak >= 3
    fnames = cell(1,2*C);
    cnt    = 1;
    for c=1:2:2*C    
        fnames{c}     = Nii.x0{cnt}(1).dat.fname;    
        fnames{c + 1} = Nii_out(cnt).dat.fname;
        cnt           = cnt + 1;
    end

    spm_check_registration(char(fnames))
end

if do_clean && (do_readwrite || (coreg && (C > 1 || numel(Nii.x{1}) > 1)))
    % Clean-up temporary files
    rmdir(dir_tmp,'s');
end
%==========================================================================
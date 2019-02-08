function GenerateTestData

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

% Where to write output
DirLowRes = './SampleData';

% Slice-profile related
DownSampling = 1/6; % Down-sampling factor, will be applied in orthogonal directions
WindowLR     = 2;
WindowHR     = 2;
Gap          = 0;

% Randomly move images rigidly?
PerturbRigid.do        = true;
PerturbRigid.transl_mx = 5;
PerturbRigid.rot_scl   = 0.05;

% Comment to use more than one observation per channel
N = 1;

% Create output directory
if  exist(DirLowRes,'dir') == 7,  rmdir(DirLowRes,'s'); end; mkdir(DirLowRes);

%--------------------------------------------------------------------------
% Get reference IXI NIfTIs
%--------------------------------------------------------------------------

dir_data = './data';
Nii_ref  = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));
C        = numel(Nii_ref);

%--------------------------------------------------------------------------
% Define projection matrices
%--------------------------------------------------------------------------

DS     = {[DownSampling 1 1; 1 DownSampling 1], ... 
          [DownSampling 1 1], ...
          [1 1 DownSampling; 1 DownSampling 1]};
  
window = {{[WindowLR WindowHR WindowHR], [WindowHR WindowLR WindowHR]}, ... 
          {[WindowLR WindowHR WindowHR]}, ...
          {[WindowHR WindowHR WindowLR], [WindowHR WindowLR WindowHR]}};

gap    = {{[Gap 0 0], [0 Gap 0]}, ... 
          {[Gap 0 0]}, ...
          {[0 0 Gap], [0 Gap 0]}};

% Sanity check  
if numel(DS) ~= C
    error('numel(DS) ~= C')
end

%--------------------------------------------------------------------------
% Simulate thick-sliced data
%--------------------------------------------------------------------------

for c=1:numel(DS) % Loop over channels
    
    % Get HR reference data
    img0 = Nii_ref(c).dat(:,:,:);
    mat0 = Nii_ref(c).mat;
    dm0  = size(img0);
        
    % Build dat object
    Nii   = {struct};
    ds    = DS{c};
    if ~exist('N','var')
        N = size(ds,1);
    end
    for n=1:N % Loop over LR images to be created
        D    = diag([ds(n,:) 1]);
        mat  = mat0/D;
        dm   = floor(D(1:3,1:3)*dm0')';
        vx   = sqrt(sum(mat(1:3,1:3).^2));   
       
        Nii{1}(n).mat     = mat;
        Nii{1}(n).dat.dim = dm;
    end
    
    dat = init_dat(Nii,mat0,dm0,window(c),gap(c));
        
    % Apply projection matrix to simulate LR data
    img = A(img0,dat);
        
    % Save LR data
    [~,nam,ext] = fileparts(Nii_ref(c).dat.fname);    
    for n=1:dat.N    
        % Save thick-sliced data        
        nfname      = fullfile(DirLowRes,['n' num2str(n) '_' nam ext]);

        % Rigidly realign the image a little bit (randomly)
        mat           = Nii{1}(n).mat;            
        dm            = Nii{1}(n).dat.dim;
        mat           = rigid_perturb(mat,dm,PerturbRigid);
        Nii{1}(n).mat = mat;
            
        % Write to NIfTI
        create_nii(nfname,img{n},Nii{1}(n).mat,[spm_type('float32') spm_platform('bigend')],'Simulated thick-sliced');
    end
end
%==========================================================================

%==========================================================================
function mat = rigid_perturb(mat,dm,PerturbRigid)
if PerturbRigid.do 
    
    B = get_rigid_basis;

    % Get random rigid perturbation
    r = zeros([1 6]);
    
    % Translation part
    mxt    = PerturbRigid.transl_mx;
    t      = randi([-mxt mxt],[3 1]) + randn([3 1]);   
    r(1:3) = t';
    
    % Rotation part
    rot_scl = PerturbRigid.rot_scl;
    r(4:6)  = rot_scl*randn(1,3); % rotation
    
    % Get matrix
    matr = spm_dexpm(r,B);    
    
    % Apply to image orientation matrix
    mat = matr*mat;
    
%     % Matrix for translating the origin
%     orig = (dm(1:3) + 1)/2;
%     off  = -orig;
%     T    = [1  0  0  off(1)
%             0  1  0  off(2)
%             0  0  1  off(3)
%             0  0  0  1     ];
% 
%     mat = T\(matr*mat)*T;
end
%==========================================================================
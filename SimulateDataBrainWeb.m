function SimulateDataBrainWeb

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

% Where to write output
DirRef = './ReferenceData/BrainWeb';
DirSim = './SimulatedData/BrainWeb';

% Randomly move images rigidly?
PerturbRigid.do        = true;
PerturbRigid.transl_mx = 1;
PerturbRigid.rot_scl   = 0;

offset = {[-2.75 1.5 -2]',[1.75 -1.5 2]',[-2 -2.5 1.5]'};
% offset = {[-1.75 1.5 -2]',[1.75 -1.5 1]',[-1 -1.5 1.5]'};

% Create output directory
DirSim3D = fullfile(DirSim,'3D');
DirSim2D = fullfile(DirSim,'2D');
if  exist(DirSim3D,'dir') ~= 7,  mkdir(DirSim3D); end
if  exist(DirSim2D,'dir') ~= 7,  mkdir(DirSim2D); end

%--------------------------------------------------------------------------
% Get noisy BrainWeb images
%--------------------------------------------------------------------------

Nii_ref = nifti(spm_select('FPList',DirRef,'^.*\pn9_rf0.nii$'));
C       = numel(Nii_ref);

%--------------------------------------------------------------------------
% Simulate data
%--------------------------------------------------------------------------

fnames = {};
for c=1:C % Loop over channels

    fname = Nii_ref(c).dat.fname;
    mat   = Nii_ref(c).mat;
    img   = Nii_ref(c).dat(:,:,:);
    dm    = size(img);
    
    fnames{end + 1} = fname;
    [~,nam,ext]     = fileparts(fname);  
    
    % Save thick-sliced data        
    nfname          = fullfile(DirSim3D,[nam ext]);
    fnames{end + 1} = nfname;

    % Rigidly realign the image a little bit (randomly)
    mat = rigid_perturb(mat,dm,PerturbRigid,offset{c});

    % Write to NIfTI
    
    % 3D
    create_nii(nfname,img,mat,[spm_type('float32') spm_platform('bigend')],'Simulated (3D)');

    % 2D
    nfname   = fullfile(DirSim2D,[nam ext]);
    mat(3,4) = 0;
    create_nii(nfname,img(:,:,floor(dm(3)/2)),mat,[spm_type('float32') spm_platform('bigend')],'Simulated (2D)');
end

spm_check_registration(char(fnames))
%==========================================================================

%==========================================================================
function mat = rigid_perturb(mat,dm,PerturbRigid,offset)
if PerturbRigid.do 
    
%     B = get_rigid_basis;
% 
%     % Get random rigid perturbation
%     r = zeros([1 6]);
%     
%     % Translation part
%     mxt    = PerturbRigid.transl_mx;
%     t      = randi([-mxt mxt],[3 1]) + randn([3 1]);   
%     r(1:3) = t';
%     
%     % Rotation part
%     rot_scl = PerturbRigid.rot_scl;
%     r(4:6)  = rot_scl*randn(1,3); % rotation
%     
%     % Get matrix
%     matr = spm_dexpm(r,B);    
%     
    % Apply to image orientation matrix
%     mat = matr*mat;
    mat(1:3,4) = mat(1:3,4) + offset;
    
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
function SimulateDataIXI

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

% Where to write output
DirRef = './ReferenceData/IXI';
DirSim = './SimulatedData/IXI';

% Slice-profile related
DownSampling = 1/6; % Down-sampling factor, will be applied in orthogonal directions
Gap          = 0;

% Comment to use more than one observation per channel
N = 1;

% Translate images a bit
offset = {[-2.75 1.5 -2]',[1.75 -1.5 2]',[-2 -2.5 1.5]'};
% offset = {[-1.75 1.5 -2]',[1.75 -1.5 1]',[-1 -1.5 1.5]'};

% Create output directory
DirSim3D = fullfile(DirSim,'3D');
DirSim2D = fullfile(DirSim,'2D');
if  exist(DirSim3D,'dir') ~= 7,  mkdir(DirSim3D); end
if  exist(DirSim2D,'dir') ~= 7,  mkdir(DirSim2D); end

%--------------------------------------------------------------------------
% Get reference IXI NIfTIs
%--------------------------------------------------------------------------

Nii_ref  = nifti(spm_select('FPList',DirRef,'^.*\.nii$'));
C        = numel(Nii_ref);

%--------------------------------------------------------------------------
% Define projection matrices
%--------------------------------------------------------------------------

DS = {[DownSampling 1 1; 1 DownSampling 1], ... 
      [DownSampling 1 1], ...
      [1 1 DownSampling; 1 DownSampling 1]};

% Sanity check  
if numel(DS) ~= C
    error('numel(DS) ~= C')
end

%--------------------------------------------------------------------------
% Simulate thick-sliced data
%--------------------------------------------------------------------------

fnames = {};
for c=1:numel(DS) % Loop over channels
    
    % Get HR reference data
    img0  = Nii_ref(c).dat(:,:,:);
    mat0  = Nii_ref(c).mat;
    dm0   = size(img0);
    
    fname           = Nii_ref(c).dat.fname;    
    fnames{end + 1} = fname;
    
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
    
    dat = init_dat('superres',Nii,mat0,dm0,[],Gap);
        
    % Apply projection matrix to simulate LR data
    img = A(img0,dat);
        
    % Save LR data
    [~,nam,ext] = fileparts(Nii_ref(c).dat.fname);    
    for n=1:dat.N    
        % Save thick-sliced data        
        nfname          = fullfile(DirSim3D,['n' num2str(n) '_' nam ext]);
        fnames{end + 1} = nfname;
        
        % Rigidly realign the image a little bit (randomly)
        mat           = Nii{1}(n).mat;          
        mat(1:3,4)    = mat(1:3,4) + offset{c};
        Nii{1}(n).mat = mat;
            
        % Write to NIfTI
        
        % 3D
        create_nii(nfname,img{n},Nii{1}(n).mat,[spm_type('float32') spm_platform('bigend')],'Simulated thick-sliced (3D)');
        
        % 2D
        nfname             = fullfile(DirSim2D,[nam ext]);
        Nii{1}(n).mat(3,4) = 0;
        create_nii(nfname,img{n}(:,:,floor(dm(3)/2)),Nii{1}(n).mat,[spm_type('float32') spm_platform('bigend')],'Simulated thick-sliced (2D)');
    end
end

spm_check_registration(char(fnames))
%==========================================================================
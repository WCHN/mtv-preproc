function SimulateDataBrainWeb

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

% Where to write output
DirRef = './ReferenceData/BrainWeb';
DirSim = './SimulatedData/BrainWeb';

ExtractSlab = true;
SlabSize    = 5;

% Translate images a bit
offset = {[-2.75 1.5 -2]',[1.75 -1.5 2]',[-2 -2.5 1]'};
% offset = {[-1.75 1.5 -2]',[1.75 -1.5 1]',[-1 -1.5 1.5]'};

% Create output directory
DirSim3D = fullfile(DirSim,'3D');
DirSim2D = fullfile(DirSim,'2D');
if  (exist(DirSim3D,'dir') == 7),  rmdir(DirSim3D,'s'); end; mkdir(DirSim3D);
if  (exist(DirSim2D,'dir') == 7),  rmdir(DirSim2D,'s'); end; mkdir(DirSim2D);

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

    if exist('offset','var')
        % Rigidly realign the image a little bit (randomly)
        mat(1:3,4) = mat(1:3,4) + offset{c};
    end
    
    % Write to NIfTI
    
    % 3D
    create_nii(nfname,img,mat,[spm_type('float32') spm_platform('bigend')],'Simulated (3D)');

    % 2D
    nfname2d = fullfile(DirSim2D,[nam ext]);
    mat(3,4) = 0;
    create_nii(nfname2d,img(:,:,floor(dm(3)/2)),mat,[spm_type('float32') spm_platform('bigend')],'Simulated (2D)');
    
    if ExtractSlab
        z      = round(dm(3)/2);
        ofname = nfname;
        bb     = [-Inf Inf; -Inf Inf; (z - SlabSize) (z + SlabSize);]';
        nfname = subvol(spm_vol(nfname),bb);
        delete(ofname);
    end
    
    fnames{end + 1} = nfname;
end

spm_check_registration(char(fnames))
%==========================================================================
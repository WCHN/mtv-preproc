function SimulateDataBrainWeb

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

% Where to write output
DirRef = './ReferenceData/BrainWeb';
DirSim = './SimulatedData/BrainWeb';

ExtractSlab3d = true;
SlabSize3d    = 5;
Plane2d       = 1;     % 1 - Sagittal, 2 - Coronal, 3 - Axial
Crop2d        = false;
CropSize2d    = 20;

% Translate images a bit
offset = {[-2.75 2.5 -1.75]',[2.25 -1.75 2.25]',[-2.25 -3.0 2.5]'};
% offset = {[-4.75 4.5 -2]',[2.75 -3.75 2]',[-2.25 -5.5 1.5]'};
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

    fname           = Nii_ref(c).dat.fname;
    fnames{end + 1} = fname;
    [~,nam,ext]     = fileparts(fname);  
    
    mat   = Nii_ref(c).mat;
    img   = Nii_ref(c).dat(:,:,:);
    dm    = size(img);    
    
    if exist('offset','var')
        % Rigidly realign the image a little bit (randomly)
        mat(1:3,4) = mat(1:3,4) + offset{c};
    end

    %---------------------
    % Write to NIfTI
    %---------------------
        
    % 3D
    nfname = fullfile(DirSim3D,[nam ext]); 
    create_nii(nfname,img,mat,[spm_type('float32') spm_platform('bigend')],'Simulated (3D)');

    % 2D
    fname2d = extract_slice(nfname,Plane2d);
    movefile(fname2d,DirSim2D);
        
    if Crop2d
        y        = round(dm(2)/2);
        ofname2d = nfname2d;
        bb       = [1 + CropSize2d dm(1) - CropSize2d; 1 + CropSize2d dm(2) - CropSize2d; -Inf Inf]';
        nfname2d = subvol(spm_vol(nfname2d),bb);
        delete(ofname2d);
    end
        
    if ExtractSlab3d
        z      = round(dm(3)/2);
        ofname = nfname;
        bb     = [-Inf Inf; -Inf Inf; (z - SlabSize3d) (z + SlabSize3d);]';
        nfname = subvol(spm_vol(nfname),bb);
        delete(ofname);
    end
    
    fnames{end + 1} = nfname;
end

spm_check_registration(char(fnames))
%==========================================================================
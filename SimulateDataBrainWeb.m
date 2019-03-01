function SimulateDataBrainWeb
% Simulate data from BrainWeb images
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
  
%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

% Where to write output
DirRef = './ReferenceData/BrainWeb';
DirSim = './SimulatedData/BrainWeb';

% If DownSampling > 0 & DownSampling < 1, creates downsampled 3D images.
% The variable NoisePrct determines the amount of noise to add to the noise
% free down-sampled BrainWeb image. Deg determines the interpolation degree.
% DownSampling = 1/3;
DownSampling = 0;
NoisePrct    = 0.15;
Deg          = 4;

% If true, reslices images to have the same size and orientation matrix.
% Also crops a little bit of the FOV.
Reslice = false;
Padding = 5;

% Translate images a bit
offset = {[3.5 -2.75 0.75]',[-3.25 0.5 -1.75]',[0.75 3.25 -3.0]'};

% Determines what plane to extract when creating 2D slices
% 1 - Sagittal, 2 - Coronal, 3 - Axial
Plane2d = 3;     

% If true, crops the 2d image so that the FOV becomes CropSize2d smalled in
% each direction
Crop2d     = false;
CropSize2d = 5;

% If true, pull out a slab from the 3D image of size 2*SlabSize3d
ExtractSlab3d = false;
SlabSize3d    = 5;

%--------------------------------------------------------------------------
% Create output directory
%--------------------------------------------------------------------------

DirSim3D = fullfile(DirSim,'3D');
DirSim2D = fullfile(DirSim,'2D');
if  (exist(DirSim3D,'dir') == 7),  rmdir(DirSim3D,'s'); end; mkdir(DirSim3D);
if  (exist(DirSim2D,'dir') == 7),  rmdir(DirSim2D,'s'); end; mkdir(DirSim2D);

%--------------------------------------------------------------------------
% Get noisy BrainWeb images
%--------------------------------------------------------------------------

if DownSampling > 0
    Nii_ref = nifti(spm_select('FPList',DirRef,'^.*\pn0_rf0.nii$'));
else
    Nii_ref = nifti(spm_select('FPList',DirRef,'^.*\pn9_rf0.nii$'));
end

C = numel(Nii_ref); % Number of channels

%--------------------------------------------------------------------------
% Simulate
%--------------------------------------------------------------------------

 % For storing filenames
fnames_in  = cell(1,C);
fnames_out = cell(1,C);

for c=1:C % Loop over channels

    fname       = Nii_ref(c).dat.fname;    
    [~,nam,ext] = fileparts(fname);  
    
    fnames_in{c} = fname;
    
    mat = Nii_ref(c).mat;
    img = Nii_ref(c).dat(:,:,:);   
    
    if DownSampling > 0
        % Down-sample image w NN interpolation
        [img,mat] = resample_img(Nii_ref(c),DownSampling,Deg);
        
        % Add noise
        msk = isfinite(img) & img > 0;
        sd  = std(img(msk));
        img = abs(img + (NoisePrct*sd)*randn(size(img)));
    end
    
    if exist('offset','var')
        % Rigidly realign the image a little bit
        mat(1:3,4) = mat(1:3,4) + offset{c};
    end    
    
    %----------------------------------------------------------------------
    % Write to NIfTI
    %----------------------------------------------------------------------
        
    % 3D
    nfname = fullfile(DirSim3D,[nam ext]); 
    create_nii(nfname,img,mat,[spm_type('float32') spm_platform('bigend')],'Simulated (3D)');

    fnames_out{c} = nfname;
    
end % End loop over channels

if Reslice
    fnames_out = reslice_imgs(fnames_out,Padding);
end

for c=1:C % Loop over channels
        
    % 2D
    fname2d     = extract_slice(fnames_out{c},Plane2d);
    [~,nam,ext] = fileparts(fname2d);  
    movefile(fname2d,DirSim2D);
    fname2d     = fullfile(DirSim2D,[nam ext]);
    
    if Crop2d
        ofname2d = fname2d;
        Nii      = nifti(ofname2d);
        dm       = Nii.dat.dim;        
        bb       = [1 + CropSize2d dm(1) - CropSize2d; 1 + CropSize2d dm(2) - CropSize2d; -Inf Inf]';
        fname2d  = subvol(spm_vol(ofname2d),bb);
        delete(ofname2d);
    end
        
    if ExtractSlab3d
        ofname = fnames_out{c};
        Nii    = nifti(ofname);
        dm     = Nii.dat.dim;
        z      = round(dm(3)/2);        
        bb     = [-Inf Inf; -Inf Inf; (z - SlabSize3d) (z + SlabSize3d);]';
        nfname = subvol(spm_vol(ofname),bb);
        delete(ofname);
        
        fnames_out{c} = nfname;
    end        
    
end % End loop

%--------------------------------------------------------------------------
% Show simulated results
%--------------------------------------------------------------------------

fnames = {};
for c=1:C
    fnames{end + 1} = fnames_in{c};
    fnames{end + 1} = fnames_out{c};
end

spm_check_registration(char(fnames))
%==========================================================================
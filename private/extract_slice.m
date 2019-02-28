function nfname = extract_slice(fname,axis_2d)
% Extract a 2D slice from a 3D volume.
% FORMAT nfname1 = create_2d_slice(fname,axis_2d,clean_up)
% fname   - Input filename
% axis_2d - Axis to extract along [axis_2d=3]
% nfname  - Filename of 2d image
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   

if nargin < 2, axis_2d  = 3; end

% Create bounding box
V  = spm_vol(fname);
dm = V.dim;
if axis_2d     == 1
    d1 = floor(dm(1)/2) + 1;
    bb = [d1 d1;-inf inf;-inf inf];   
elseif axis_2d == 2
    d1 = floor(dm(2)/2) + 1;
    bb = [-inf inf;d1 d1;-inf inf];
elseif axis_2d == 3 
    d1 = floor(dm(3)/2) + 1;
    bb = [-inf inf;-inf inf;d1 d1];
end                

% Crop according to bounding-box
nfname = subvol(V,bb');

% Make sure 1D plane is in z dimension
Nii  = nifti(nfname);
mat  = Nii.mat;
img  = Nii.dat(:,:,:);

if axis_2d == 1 || axis_2d == 2
    % Permute image data and apply permutation matrix to orientation matrix
    if axis_2d == 1
        img = permute(img,[2 3 1]);            
        P   = [0 1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];
    else
        img = permute(img,[1 3 2]);        
        P   = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
    end   

    mat     = P*mat*P';    
end

% Set zero translation in z-direction
mat(3,4) = 0;

% Overwrite image data
VO             = spm_vol(nfname);
dm             = [size(img) 1];
VO.dim(1:3)    = dm(1:3);        
VO.mat         = mat;
VO             = spm_create_vol(VO);        
Nii            = nifti(VO.fname);    
Nii.dat(:,:,:) = img; 
%==========================================================================
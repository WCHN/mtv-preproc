% clear; clc;

% Down-sampling factor, will be applied in orthogonal directions
DownSampling = 1/6;

% Get NIfTIs
dir_data = './data';
Nii_ref  = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));
C        = numel(Nii_ref);

% Simulate thick-sliced data
ds     = [DownSampling 1 1; 1 DownSampling 1; 1 1 DownSampling];
Nii_in = nifti;
for c=1:C
    img  = Nii_ref(c).dat(:,:,:);
    mat0 = Nii_ref(c).mat;
    dm0  = size(img);
    D    = diag([ds(c,:) 1]);
    mat  = mat0/D;
    dm   = floor(D(1:3,1:3)*dm0')';
    
    Nii_dat.mat     = mat;
    Nii_dat.dat.dim = dm;
    dat             = init_dat(Nii_dat,mat0,dm0);
    
    % Apply projection matrix
    img = A(img,dat);
    
    % Rescale intensities
    vx0 = sqrt(sum(mat0(1:3,1:3).^2)); 
    vx  = sqrt(sum(mat(1:3,1:3).^2)); 
    scl = prod(vx0./vx);
    img = scl*img{1};
        
    % Save thick-sliced data
    [~,nam,ext] = fileparts(Nii_ref(c).dat.fname);
    nfname      = fullfile(dir_data,['ds_' nam ext]);
    
    VO          = spm_vol(Nii_ref(c).dat.fname);
    VO.fname    = nfname;
    VO.dim(1:3) = dm(1:3);    
    VO.mat      = mat;
    VO          = spm_create_vol(VO);
        
    Nii               = nifti(VO.fname);
    nii.dat.scl_slope = max(img(:))/1600;
    Nii.dat(:,:,:)    = img;
end
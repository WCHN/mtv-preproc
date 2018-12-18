function GenerateTestData

% Down-sampling factor, will be applied in orthogonal directions
DownSampling = 1/6;
DirLowRes    = 'LowResData';

% Create output directory
if  exist(DirLowRes,'dir') == 7,  rmdir(DirLowRes,'s'); end; mkdir(DirLowRes);

% Get reference IXI NIfTIs
dir_data = './data';
Nii_ref  = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));
C        = numel(Nii_ref);

% Set downsampling
DS = {[DownSampling 1 1; 1 DownSampling 1], ... 
      [1 DownSampling 1], ...
      [1 1 DownSampling; 1 DownSampling 1]};
  
% Sanity check  
if numel(DS) ~= C
    error('numel(DS) ~= C')
end

% Simulate thick-sliced data
for c=1:numel(DS) % Loop over channels
    
    % Get HR reference data
    img0 = Nii_ref(c).dat(:,:,:);
    mat0 = Nii_ref(c).mat;
    dm0  = size(img0);
        
    % Build dat object
    Nii = {struct};
    ds  = DS{c};
    N   = size(ds,1);
    for n=1:N % Loop over LR images to be created
        D    = diag([ds(n,:) 1]);
        mat  = mat0/D;
        dm   = floor(D(1:3,1:3)*dm0')';

        Nii{1}(n).mat     = mat;
        Nii{1}(n).dat.dim = dm;
    end
    
    dat = init_dat(Nii,mat0,dm0);
        
    % Apply projection matrix to simulate LR data
    img = A(img0,dat);

    % Save LR data
    [~,nam,ext] = fileparts(Nii_ref(c).dat.fname);    
    for n=1:dat.N    
        % Save thick-sliced data        
        nfname      = fullfile(DirLowRes,['ds_n' num2str(n) '_' nam ext]);

        % Write to NIfTI
        create_nii(nfname,img{n},Nii{1}(n).mat,[spm_type('float32') spm_platform('bigend')],'Simulated thick-sliced');
    end
end
clear;

dir_data = './ReferenceData/BrainWeb';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Split noisy and references
Nii_ref   = Nii([1 3 5]);
Nii_noisy = Nii([2 4 6]);
C         = 3;

% Parameters
WorkersParfor   = Inf;
Verbose         = 2;
IterMax         = 20;
ADMMStepSize    = 0;
CoRegister      = false;
EstimateRigid   = false;

% Grid-searches if vector
RegScaleDenoisingMRI = 1:0.5:10; 

psnrs = zeros([1 numel(RegScaleDenoisingMRI) + 1]);
ssims = zeros([1 numel(RegScaleDenoisingMRI) + 1]);

ref  = cat(3,single(Nii_ref(1).dat(:,:,:)),single(Nii_ref(2).dat(:,:,:)),single(Nii_ref(3).dat(:,:,:)));
nois = cat(3,single(Nii_noisy(1).dat(:,:,:)),single(Nii_noisy(2).dat(:,:,:)),single(Nii_noisy(3).dat(:,:,:)));

psnrs(1) = get_psnr(nois(:),ref(:));
ssims(1) = ssim(nois,ref);
clear ref nois

for r=1:numel(RegScaleDenoisingMRI)
    
    OutputDirectory = ['./Output/DenoiseBrainWeb/' num2str(RegScaleDenoisingMRI(r))];
    
    % Denoise
    Nii_den = spm_mtv_preproc('InputImages',Nii_noisy,'Verbose',Verbose, ...
                              'WorkersParfor',WorkersParfor, ...
                              'RegScaleDenoisingMRI',RegScaleDenoisingMRI(r),'IterMax',IterMax,'Reference',Nii_ref, ...
                              'ADMMStepSize',ADMMStepSize,'OutputDirectory',OutputDirectory, ...
                              'CoRegister',CoRegister,'EstimateRigid',EstimateRigid);
    
    % Compute metrics
    ref  = cat(3,single(Nii_ref(1).dat(:,:,:)),single(Nii_ref(2).dat(:,:,:)),single(Nii_ref(3).dat(:,:,:)));
    den = cat(3,single(Nii_den(1).dat(:,:,:)),single(Nii_den(2).dat(:,:,:)),single(Nii_den(3).dat(:,:,:)));

    psnrs(1 + r) = get_psnr(den(:),ref(:));
    ssims(1 + r) = ssim(den,ref);
    clear ref den                          
end                      
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
RegScaleDenoisingMRI = 1:0.5:20; 

[psnrs(1),ssims(1)] = compute_image_metrics(Nii_noisy,Nii_ref);

for r=1:numel(RegScaleDenoisingMRI)
    
    OutputDirectory = ['./Output/DenoiseBrainWeb/' num2str(RegScaleDenoisingMRI(r))];
    
    % Denoise
    Nii_den = spm_mtv_preproc('InputImages',Nii_noisy,'Verbose',Verbose, ...
                              'WorkersParfor',WorkersParfor, ...
                              'RegScaleDenoisingMRI',RegScaleDenoisingMRI(r),'IterMax',IterMax,'Reference',Nii_ref, ...
                              'ADMMStepSize',ADMMStepSize,'OutputDirectory',OutputDirectory, ...
                              'CoRegister',CoRegister,'EstimateRigid',EstimateRigid);
    
    % Compute metrics
    [psnrs(1 + r),ssims(1 + r)] = compute_image_metrics(Nii_den,Nii_ref);
end                      

%%
figure(666);
subplot(121)
plot([0 RegScaleDenoisingMRI], psnrs,'b-'); title('PSNR')
subplot(122)
plot([0 RegScaleDenoisingMRI], ssims,'r-'); title('SSIM')
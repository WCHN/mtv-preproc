clear;

dir_data = './ReferenceData/BrainWeb';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Split noisy and references
Nii_ref   = Nii([1 3 5]);
Nii_noisy = Nii([2 4 6]);
C         = 3;

% Parameters
WorkersParfor = Inf;
Verbose       = 0;
IterMax       = 20;
ADMMStepSize  = 0;
CoRegister    = false;
EstimateRigid = false;
Reference     = [];%Nii_ref;

%% Run multi-channel denoising
RegScaleDenoisingMRI = 1:0.25:15; % Grid-searches if vector [5]

[psnrs(1),ssims(1)] = compute_image_metrics(Nii_noisy,Nii_ref);

for r=1:numel(RegScaleDenoisingMRI)
    fprintf('%g ',RegScaleDenoisingMRI(r))
    
    OutputDirectory = ['./Output/DenoiseBrainWeb/' num2str(RegScaleDenoisingMRI(r))];
    
    % Denoise
    tic;
    Nii_den = spm_mtv_preproc('InputImages',Nii_noisy,'Verbose',Verbose, ...
                              'WorkersParfor',WorkersParfor, ...
                              'RegScaleDenoisingMRI',RegScaleDenoisingMRI(r),'IterMax',IterMax,'Reference',Reference, ...
                              'ADMMStepSize',ADMMStepSize,'OutputDirectory',OutputDirectory, ...
                              'CoRegister',CoRegister,'EstimateRigid',EstimateRigid);
    t(r) = toc;
    
    % Compute metrics
    [psnrs(1 + r),ssims(1 + r)] = compute_image_metrics(Nii_den,Nii_ref);
end                      
fprintf('\n')

% Plot results
figure(666);
subplot(221)
plot([0 RegScaleDenoisingMRI], psnrs,'b-'); title('PSNR')
subplot(222)
plot([0 RegScaleDenoisingMRI], ssims,'r-'); title('SSIM')
subplot(223)
plot(RegScaleDenoisingMRI, t,'g-'); title('Runtime')

%% Compare to other method(s)
% https://sites.google.com/site/pierrickcoupe/softwares/denoising-for-medical-imaging/mri-denoising/mri-denoising-software
% ONLM filter




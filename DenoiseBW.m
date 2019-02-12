clear;

dir_data = './SimulatedData/BrainWeb';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Parameters
WorkersParfor        = 0;
Verbose              = 2;
IterMax              = 20;
ADMMStepSize         = 0;
CoRegister           = false;
EstimateRigid        = true;
RegScaleDenoisingMRI = 5;
OutputDirectory      = 'Output/DenoiseRigid';
DecreasingReg        = true;
IterImage            = 2;

Nii_den = spm_mtv_preproc('InputImages',Nii,'Verbose',Verbose, ...
                          'WorkersParfor',WorkersParfor, ...
                          'RegScaleDenoisingMRI',RegScaleDenoisingMRI,'IterMax',IterMax, ...
                          'ADMMStepSize',ADMMStepSize,'OutputDirectory',OutputDirectory, ...
                          'CoRegister',CoRegister,'EstimateRigid',EstimateRigid, ...
                          'DecreasingReg',DecreasingReg,'IterImage',IterImage);

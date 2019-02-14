clear;

dir_data = './SimulatedData/BrainWeb/2D';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Parameters
WorkersParfor   = Inf;
Verbose         = 3;
OutputDirectory = 'Output/DenoiseRigid';

CoRegister           = false;
EstimateRigid        = true;
ADMMStepSize         = 0;
IterMax              = 30;
Tolerance            = 0;
MeanCorrectRigid     = true;
RegScaleDenoisingMRI = 6;
IterImage            = 3;

Nii_den = spm_mtv_preproc('InputImages',Nii,'Verbose',Verbose, ...
                          'WorkersParfor',WorkersParfor, ...                          
                          'OutputDirectory',OutputDirectory, ...
                          'CoRegister',CoRegister,'EstimateRigid',EstimateRigid, ...
                          'ADMMStepSize',ADMMStepSize,'MeanCorrectRigid',MeanCorrectRigid, ...
                          'IterMax',IterMax,'Tolerance',Tolerance,'IterImage',IterImage, ...
                          'RegScaleDenoisingMRI',RegScaleDenoisingMRI);

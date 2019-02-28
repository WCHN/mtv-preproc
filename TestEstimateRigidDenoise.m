clear;

dir_data = './SimulatedData/BrainWeb/2D';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Parameters
WorkersParfor   = Inf;
Verbose         = 3;
OutputDirectory = 'Output/DenoiseRigid';

CoRegister           = false;

DecreasingReg        = true;
EstimateRigid        = true;
PaddingBB            = 5;

ADMMStepSize         = 0;
IterMax              = 30;
Tolerance            = 1e-4;
MeanCorrectRigid     = true;
RegScaleDenoisingMRI = 10;

IterImage            = 12;
IterGaussNewtonRigid = 3;
IterGaussNewtonImage = 1;

Nii_den = spm_mtv_preproc('InputImages',Nii,'Verbose',Verbose, ...
                          'WorkersParfor',WorkersParfor, ...                          
                          'OutputDirectory',OutputDirectory, ...
                          'CoRegister',CoRegister,'EstimateRigid',EstimateRigid, ...
                          'ADMMStepSize',ADMMStepSize,'MeanCorrectRigid',MeanCorrectRigid, ...
                          'IterMax',IterMax,'Tolerance',Tolerance,'IterImage',IterImage, ...
                          'RegScaleDenoisingMRI',RegScaleDenoisingMRI, ...
                          'IterGaussNewtonRigid',IterGaussNewtonRigid, ...
                          'IterGaussNewtonImage',IterGaussNewtonImage, ...
                          'DecreasingReg',DecreasingReg,'PaddingBB',PaddingBB);

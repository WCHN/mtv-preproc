clear;

%--------------------------------------------------------------------------
% Read input data
%--------------------------------------------------------------------------

dir_data = './SimulatedData/BrainWeb/ThickSliced/2D';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

%--------------------------------------------------------------------------
% Algorithm parameters
%--------------------------------------------------------------------------

WorkersParfor   = Inf;
Verbose         = 3;  % 0 = silent, 1 = text, 2 = +figures, 3 = +spm_chec_reg
OutputDirectory = 'Output/TestEstimateRigidSuperres';

CoRegister           = true;
DecreasingReg        = true;
EstimateRigid        = true;

VoxelSize            = 1;
PaddingBB            = 0;
ADMMStepSize         = 0;
IterMax              = 100;
Tolerance            = 1e-4;
MeanCorrectRigid     = true;
RegScaleDenoisingMRI = 40;

IterImage            = 1;
IterGaussNewtonRigid = 3;
IterGaussNewtonImage = 1;

%--------------------------------------------------------------------------
% Run algorithm
%--------------------------------------------------------------------------

Nii_den = spm_mtv_preproc('InputImages',Nii,'Verbose',Verbose,'Method','superres', ...
                          'WorkersParfor',WorkersParfor, ...                          
                          'OutputDirectory',OutputDirectory, ...
                          'CoRegister',CoRegister,'EstimateRigid',EstimateRigid, ...
                          'ADMMStepSize',ADMMStepSize,'MeanCorrectRigid',MeanCorrectRigid, ...
                          'IterMax',IterMax,'Tolerance',Tolerance,'IterImage',IterImage, ...
                          'RegScaleDenoisingMRI',RegScaleDenoisingMRI, ...
                          'IterGaussNewtonRigid',IterGaussNewtonRigid, ...
                          'IterGaussNewtonImage',IterGaussNewtonImage, ...
                          'DecreasingReg',DecreasingReg,'PaddingBB',PaddingBB, ...
                          'VoxelSize',VoxelSize);            
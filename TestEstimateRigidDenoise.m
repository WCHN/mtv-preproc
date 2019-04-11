clear;

%--------------------------------------------------------------------------
% Read input data
%--------------------------------------------------------------------------

NoisePrct = 0;
Offset    = {[5.5; 4.25; 0],[0; 0; 0],[-4.5; -3.75; 0]};
Rotation  = {[0; 0; 0],[0; 0; 0],[0; 0; 0]};

[Nii_sim3d,Nii_sim2d] = SimulateDataBrainWeb(NoisePrct,Offset,Rotation);
Nii                   = Nii_sim2d;

%--------------------------------------------------------------------------
% Algorithm parameters
%--------------------------------------------------------------------------

WorkersParfor   = Inf;
Verbose         = 2;  % 0 = silent, 1 = text, 2 = +figures, 3 = +spm_chec_reg
OutputDirectory = 'Output/TestEstimateRigidDenoise';
ReadWrite       = false;

CoRegister           = false;
DecreasingReg        = true;
EstimateRigid        = true;

VoxelSize            = sqrt(sum(Nii(1).mat(1:3,1:3).^2));
PaddingBB            = 0;
ADMMStepSize         = 0;
IterMax              = 100;
Tolerance            = 1e-3;
MeanCorrectRigid     = true;
RegScaleDenoisingMRI = 30;

IterImage            = 5;
IterGaussNewtonRigid = 5;
IterGaussNewtonImage = 1;

%--------------------------------------------------------------------------
% Run algorithm
%--------------------------------------------------------------------------

Nii_den = spm_mtv_preproc('InputImages',Nii,'Verbose',Verbose, ...
                          'WorkersParfor',WorkersParfor, ...                          
                          'OutputDirectory',OutputDirectory, ...
                          'CoRegister',CoRegister,'EstimateRigid',EstimateRigid, ...
                          'ADMMStepSize',ADMMStepSize,'MeanCorrectRigid',MeanCorrectRigid, ...
                          'IterMax',IterMax,'Tolerance',Tolerance,'IterImage',IterImage, ...
                          'RegScaleDenoisingMRI',RegScaleDenoisingMRI, ...
                          'IterGaussNewtonRigid',IterGaussNewtonRigid, ...
                          'IterGaussNewtonImage',IterGaussNewtonImage, ...
                          'DecreasingReg',DecreasingReg,'PaddingBB',PaddingBB, ...
                          'VoxelSize',VoxelSize,'ReadWrite',ReadWrite);
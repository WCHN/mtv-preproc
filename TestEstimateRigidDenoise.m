clear;

dir_data = './SimulatedData/BrainWeb';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Parameters
WorkersParfor   = Inf;
Verbose         = 2;
OutputDirectory = 'Output/DenoiseRigid';

CoRegister    = true;
EstimateRigid = true;

Nii_den = spm_mtv_preproc('InputImages',Nii,'Verbose',Verbose, ...
                          'WorkersParfor',WorkersParfor, ...                          
                          'OutputDirectory',OutputDirectory, ...
                          'CoRegister',CoRegister,'EstimateRigid',EstimateRigid);

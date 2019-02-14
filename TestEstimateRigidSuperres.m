clear;

dir_data = './SimulatedData/IXI/2D';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));


Verbose         = 2;
WorkersParfor   = Inf;
OutputDirectory = 'Output/SuperresRigid';

IterMax             = 30;
VoxelSize           = 1;
ADMMStepSize        = 0;
RegScaleSuperResMRI = 6;
DecreasingReg       = true;
CoRegister          = false;
IterImage           = 3;
EstimateRigid       = false;

spm_mtv_preproc('InputImages',Nii,'Method','superres', ...
                'VoxelSize',VoxelSize,'Verbose',Verbose, ...
                'RegScaleSuperResMRI',RegScaleSuperResMRI, ...                
                'ADMMStepSize',ADMMStepSize,'WorkersParfor',WorkersParfor, ...
                'IterMax',IterMax, ...
                'DecreasingReg',DecreasingReg,'CoRegister',CoRegister, ...
                'EstimateRigid',EstimateRigid,'IterImage',IterImage, ...
                'OutputDirectory',OutputDirectory); 
clear;

dir_data = './SimulatedData/IXI/3D';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));


Verbose         = 3;
WorkersParfor   = Inf;
OutputDirectory = 'Output/SuperresRigid';

IterMax             = 40;
VoxelSize           = 1;
ADMMStepSize        = 0;
RegScaleSuperResMRI = 6;
DecreasingReg       = true;
CoRegister          = false;
IterImage           = 1;
EstimateRigid       = true;

spm_mtv_preproc('InputImages',Nii,'Method','superres', ...
                'VoxelSize',VoxelSize,'Verbose',Verbose, ...
                'RegScaleSuperResMRI',RegScaleSuperResMRI, ...                
                'ADMMStepSize',ADMMStepSize,'WorkersParfor',WorkersParfor, ...
                'IterMax',IterMax, ...
                'DecreasingReg',DecreasingReg,'CoRegister',CoRegister, ...
                'EstimateRigid',EstimateRigid,'IterImage',IterImage, ...
                'OutputDirectory',OutputDirectory); 
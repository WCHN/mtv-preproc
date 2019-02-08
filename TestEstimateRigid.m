clear;

dir_data = './SampleData';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

VoxelSize           = 1;
Verbose             = 3;
WorkersParfor       = Inf;

IterMax             = 20;
Tolerance           = 0;

ADMMStepSize        = 0;
RegScaleSuperResMRI = 20;
DecreasingReg       = true;
CoRegister          = false;
EstimateRigid       = true;
IterImage           = 3;

spm_mtv_preproc('InputImages',Nii,'Method','superres', ...
                'VoxelSize',VoxelSize,'Verbose',Verbose, ...
                'RegScaleSuperResMRI',RegScaleSuperResMRI, ...                
                'ADMMStepSize',ADMMStepSize,'WorkersParfor',WorkersParfor, ...
                'Tolerance',Tolerance,'IterMax',IterMax, ...
                'DecreasingReg',DecreasingReg,'CoRegister',CoRegister, ...
                'EstimateRigid',EstimateRigid,'IterImage',IterImage);
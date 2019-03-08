clear;

%--------------------------------------------------------------------------
% Read input data
%--------------------------------------------------------------------------

dir_data = './SimulatedData/BrainWeb/3D';
Nii_sim  = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

Nii          = struct;
[Nii,C,is3d] = parse_input_data(Nii,Nii_sim,'denoise');

dir_tmp  = 'Temp/TestMIReg';
if  exist(dir_tmp,'dir') == 7  
    rmdir(dir_tmp,'s'); 
end
mkdir(dir_tmp); 

Nii.x = copy_ims(Nii.x,dir_tmp);
Nii.x = coreg_ims(Nii.x);

fnames = cell(1,2*C);
cnt    = 1;
for c=1:2:2*C    
    fnames{c}     = Nii_sim(cnt).dat.fname;    
    fnames{c + 1} = Nii.x{cnt}.dat.fname;
    cnt           = cnt + 1;
end

spm_check_registration(char(fnames))
function rigidly_realign(fname,offset,rotation)
% Rigidly (translation and rotation) modify an orientation matrix
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

Nii  = nifti(fname);
dm   = Nii.dat.dim;
is2d = numel(dm) == 2;

if is2d
    q = [offset(1) offset(2) 0 0 0 rotation(3)];
else
    q = [offset(1) offset(2) offset(3) rotation(1) rotation(2) rotation(3)];
end
M  = spm_matrix(q(:)');
MM = spm_get_space(fname);
spm_get_space(fname, M\MM);
%==========================================================================
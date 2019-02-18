function Nii = alloc_aux_vars(Nii,do_readwrite,C,dm,mat,dir_tmp)
% Allocate MTV model auxiliary variables
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if do_readwrite
    % Read/write temporary variables from disk (stored as NIfTIs)
    Nii.y = nifti;
    Nii.u = nifti;
    Nii.w = nifti;
    Nii.H = nifti;
else
    % Keep temporary variables in memory
    Nii.y = struct;
    Nii.u = struct;
    Nii.w = struct;
    Nii.H = struct;
end

for c=1:C
    if do_readwrite
        fname_y = fullfile(dir_tmp,['y' num2str(c) '.nii']);
        fname_u = fullfile(dir_tmp,['u' num2str(c) '.nii']); 
        fname_w = fullfile(dir_tmp,['w' num2str(c) '.nii']);
        fname_H = fullfile(dir_tmp,['H' num2str(c) '.nii']);

        create_nii(fname_y,zeros(dm,      'single'),mat,[spm_type('float32') spm_platform('bigend')],'y');
        create_nii(fname_u,zeros([dm 3 2],'single'),mat,[spm_type('float32') spm_platform('bigend')],'u');
        create_nii(fname_w,zeros([dm 3 2],'single'),mat,[spm_type('float32') spm_platform('bigend')],'w');
        create_nii(fname_H,zeros(dm,      'single'),mat,[spm_type('float32') spm_platform('bigend')],'H');

        Nii.y(c) = nifti(fname_y);
        Nii.u(c) = nifti(fname_u);
        Nii.w(c) = nifti(fname_w);
        Nii.H(c) = nifti(fname_H);
    else
        Nii.y(c).dat = zeros(dm,      'single');
        Nii.u(c).dat = zeros([dm 3 2],'single');
        Nii.w(c).dat = zeros([dm 3 2],'single');
        Nii.H(c).dat = zeros(dm,      'single');
    end
end
%==========================================================================
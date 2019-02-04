function [Nii_y,Nii_u,Nii_w] = alloc_aux_vars(do_readwrite,C,dm,mat,dir_tmp)
% Allocate MTV model auxiliary variables
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if do_readwrite
    % Read/write temporary variables from disk (stored as NIfTIs)
    Nii_y = nifti;
    Nii_u = nifti;
    Nii_w = nifti;
else
    % Keep temporary variables in memory
    Nii_y = struct;
    Nii_u = struct;
    Nii_w = struct;
end
for c=1:C
    if do_readwrite
        fname_y = fullfile(dir_tmp,['y' num2str(c) '.nii']);
        fname_u = fullfile(dir_tmp,['u' num2str(c) '.nii']); 
        fname_w = fullfile(dir_tmp,['w' num2str(c) '.nii']);

        create_nii(fname_y,zeros(dm,    'single'),mat,[spm_type('float32') spm_platform('bigend')],'y');
        create_nii(fname_u,zeros([dm 3 2],'single'),mat,[spm_type('float32') spm_platform('bigend')],'u');
        create_nii(fname_w,zeros([dm 3 2],'single'),mat,[spm_type('float32') spm_platform('bigend')],'w');

        Nii_y(c) = nifti(fname_y);
        Nii_u(c) = nifti(fname_u);
        Nii_w(c) = nifti(fname_w);
    else
        Nii_y(c).dat = zeros(dm,    'single');
        Nii_u(c).dat = zeros([dm 3 2],'single');
        Nii_w(c).dat = zeros([dm 3 2],'single');
    end
end
%==========================================================================
function Nii = alloc_aux_vars(Nii,do_readwrite,dm,mat,use_projmat,p)
% Allocate MTV model auxiliary variables
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C             = numel(Nii.x);
dir_tmp       = p.Results.TemporaryDirectory;
EstimateBias  = p.Results.EstimateBias;

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
Nii.b = cell(1,C);

for c=1:C
    if do_readwrite
        fname_y = fullfile(dir_tmp,['y' num2str(c) '.nii']);
        fname_u = fullfile(dir_tmp,['u' num2str(c) '.nii']); 
        fname_w = fullfile(dir_tmp,['w' num2str(c) '.nii']);
        if use_projmat
            fname_H = fullfile(dir_tmp,['H' num2str(c) '.nii']);
        end
        
        create_nii(fname_y,zeros(dm,      'single'),mat,[spm_type('float32') spm_platform('bigend')],'y');
        create_nii(fname_u,zeros([dm 3 2],'single'),mat,[spm_type('float32') spm_platform('bigend')],'u');
        create_nii(fname_w,zeros([dm 3 2],'single'),mat,[spm_type('float32') spm_platform('bigend')],'w');
        if use_projmat
            create_nii(fname_H,zeros(dm,      'single'),mat,[spm_type('float32') spm_platform('bigend')],'H');
        end
        
        Nii.y(c) = nifti(fname_y);
        Nii.u(c) = nifti(fname_u);
        Nii.w(c) = nifti(fname_w);
        
        if use_projmat
            Nii.H(c) = nifti(fname_H);
        else
            Nii.H(c) = nifti;
        end        
        
        Nii.b{c} = nifti;
        for n=1:numel(Nii.x{c})
            if EstimateBias    
                fname_b     = fullfile(dir_tmp,['bf' num2str(c) num2str(n) '.nii']);  
                create_nii(fname_b,zeros(dm,'single'),mat,[spm_type('float32') spm_platform('bigend')],'bf');
                Nii.b{c}(n) = nifti(fname_b);
            else
                Nii.b{c}(n).dat = nifti;
            end
        end
    else
        Nii.y(c).dat = zeros(dm,      'single');
        Nii.u(c).dat = zeros([dm 3 2],'single');
        Nii.w(c).dat = zeros([dm 3 2],'single');
        if use_projmat
            Nii.H(c).dat = zeros(dm,      'single');
        else
            Nii.H(c).dat = 0;
        end
        
        for n=1:numel(Nii.x{c})
            if EstimateBias                    
                Nii.b{c}(n).dat = zeros(dm,'single');
            else
                Nii.b{c}(n).dat = 0;
            end
        end
    end
end
%==========================================================================
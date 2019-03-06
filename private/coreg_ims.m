function Nii = coreg_ims(Nii)
% Alternately co-register images so that each channel acts as the reference
% image
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(Nii); % Number of image channels

if C == 1 && numel(Nii{1}) == 1
    % Only one image, no need to co-register
    return;
end

% Collect all image in spm_vol object (V)
V   = spm_vol;
cnt = 1;
for c=1:C
    N = numel(Nii{c});
    for n=1:N
        f      = Nii{c}(n).dat.fname;
        V(cnt) = spm_vol(f);
        cnt    = cnt + 1;
    end
end

% Set registration options
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];

% Start co-registration
for ref=1:1%numel(V) % Iterate over all images
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {V(ref).fname}; % Reference image to register to
    
    source = 1:numel(V);
    source = source(source ~= ref);
    for i=source % Iterate over all images, except reference
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {V(i).fname}; 

        spm_jobman('run',matlabbatch);

        V(i) = spm_vol(V(i).fname);
    end
end

% Update NIfTIs, to make sure they have up-to-date orientation matrices
for c=1:C
    N = numel(Nii{c});
    for n=1:N        
        Nii{c}(n) = nifti(Nii{c}(n).dat.fname);
    end
end
%==========================================================================

%==========================================================================
function [c,n] = get_ix(Nii,i)
C   = numel(Nii);
cnt = 1;
for c=1:C
    N = numel(Nii{c});
    for n=1:N
        if i == cnt, return; end
            
        cnt = cnt + 1;
    end
end
%==========================================================================
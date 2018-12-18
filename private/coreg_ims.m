function Nii = coreg_ims(Nii,dir_write)
% Co-register images
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(Nii);
if C == 1 && numel(Nii{1}) == 1
    return;
end

V   = spm_vol;
cnt = 1;
for c=1:C
    N = numel(Nii{c});
    for n=1:N
        f           = Nii{c}(n).dat.fname;
        [~,nam,ext] = fileparts(f);
        nf          = fullfile(dir_write,[nam ext]);

        copyfile(f,nf);

        Nii{c}(n) = nifti(nf);
        V(cnt)    = spm_vol(nf);
        cnt       = cnt + 1;
    end
end

% Get image with smallest voxel size and pick this image as reference
prod_vx = zeros(1,C);
for i=1:numel(V)
    vx         = sqrt(sum(V(i).mat(1:3,1:3).^2));    
    prod_vx(i) = prod(vx);
end
[~,ref] = min(prod_vx);

% Set options
matlabbatch{1}.spm.spatial.coreg.estimate.ref               = {V(ref).fname};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];

% Co-register
source = 1:numel(V);
source = source(source ~= ref);
for i=source
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {V(i).fname}; 

    spm_jobman('run',matlabbatch);
    
    [c,n]     = get_ix(Nii,i);
    Nii{c}(n) = nifti(V(i).fname);
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
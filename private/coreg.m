function Nii = coreg(Nii,dir_tmp)
% Co-register images
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(Nii);
if C == 1
    return;
end

V = spm_vol;
for c=1:C
    f           = Nii(c).dat.fname;
    [~,nam,ext] = fileparts(f);
    nf          = fullfile(dir_tmp,[nam ext]);
    
    copyfile(f,nf);
    
    Nii(c) = nifti(nf);
    V(c)   = spm_vol(nf);
end


% Get image with smallest voxel size and pick this image as reference
prod_vx = zeros(1,C);
for c=1:C
    vx         = sqrt(sum(V(c).mat(1:3,1:3).^2));    
    prod_vx(c) = prod(vx);
end
[~,ref] = min(prod_vx);

% Set options
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {V(ref).fname};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% Co-register
source = 1:C;
source = source(source ~= ref);
for c=source
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {V(c).fname}; 

    spm_jobman('run',matlabbatch);
    
    Nii(c) = nifti(V(c).fname);
end
%==========================================================================
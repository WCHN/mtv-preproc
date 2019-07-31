function Nii = copy_ims(Nii)
% Copy Niis to dir_write and update nifti struct.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(Nii);

for c=1:C
    N = numel(Nii{c});
    for n=1:N
        f             = Nii{c}(n).dat.fname;
        [pth,nam,ext] = fileparts(f);
        nf            = fullfile(pth,['cpy' nam ext]);

        copyfile(f,nf);

        Nii{c}(n) = nifti(nf);
    end
end
%==========================================================================
function fnames = reslice_imgs(fnames,padding)
% Reslices input images to have the same size and orientation matrix.
% Also crops a little bit of the FOV.
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 2, padding = 0; end

% Get bounding-box of SPM template
PthTemplateSPM = fullfile(spm('dir'),'tpm','TPM.nii,');
V              = spm_vol(PthTemplateSPM);
bb             = world_bb(V(1));

% Pad BB
bb(1,:) = bb(1,:) - padding;
bb(2,:) = bb(2,:) + padding;

% Collect spm_vol
C = numel(fnames);
V = spm_vol;
for c=1:C
    V(c) = spm_vol(fnames{c});
end

% Reslice
nV = reslice_to_bb(V,bb,[],'res_',0);

for c=1:C
    fnames{c} = nV(c).fname;
    delete(V(c).fname);
end
%==========================================================================

%==========================================================================
function Vo = reslice_to_bb(Vi,bb,vx,prefix,deg)

% reslice images one-by-one
Vo = spm_vol;
c  = 1;
for V=Vi
    % (copy to allow defaulting of NaNs differently for each volume)
    if isempty(vx)
        vx = sqrt(sum(V.mat(1:3,1:3).^2));        
    end
    
    voxdim = vx;
    % default voxdim to current volume's voxdim, (from mat parameters)
    if any(isnan(voxdim))
        vprm = spm_imatrix(V.mat);
        vvoxdim = vprm(7:9);
        voxdim(isnan(voxdim)) = vvoxdim(isnan(voxdim));
    end
    voxdim = voxdim(:)';

    mn = bb(1,:);
    mx = bb(2,:);
    % default BB to current volume's
    if any(isnan(bb(:)))
        vbb = world_bb(V);
        vmn = vbb(1,:);
        vmx = vbb(2,:);
        mn(isnan(mn)) = vmn(isnan(mn));
        mx(isnan(mx)) = vmx(isnan(mx));
    end

    if sum(bb(:,3)) == 0
        offset = 20;
        mn(2)  =  mn(2) + offset;
        mx(2)  =  mx(2) + offset;
    end
    
    % voxel [1 1 1] of output should map to BB mn
    % (the combination of matrices below first maps [1 1 1] to [0 0 0])
    mat = spm_matrix([mn 0 0 0 voxdim])*spm_matrix([-1 -1 -1]);
    % voxel-coords of BB mx gives number of voxels required
    % (round up if more than a tenth of a voxel over)
    imgdim = ceil(mat \ [mx 1]' - 0.1)';

    % output image
    Vo(c)       = V;
    [pth,nam,ext] = fileparts(V.fname);
    Vo(c).fname      = fullfile(pth,[prefix nam ext]);
    Vo(c).dim(1:3)   = imgdim(1:3);
    Vo(c).mat        = mat;
    Vo(c) = spm_create_vol(Vo(c));
    for i = 1:imgdim(3)
        
        D = diag([-1 1 1 1]);
        if det(V.mat(1:3,1:3)) < 0
            D = diag([1 1 1 1]);
        end

        M = inv(spm_matrix([0 0 -i])*inv(Vo(c).mat)*(D*V.mat));
        img = spm_slice_vol(V, M, imgdim(1:2), deg);
        
        spm_write_plane(Vo(c), img, i);
    end
    
    c = c + 1;
end
%==========================================================================

%==========================================================================
function bb = world_bb(V)
%  world-bb -- get bounding box in world (mm) coordinates

d = V.dim(1:3);
% corners in voxel-space
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
% corners in world-space
tc = V.mat(1:3,1:4)*c;

% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];
%==========================================================================
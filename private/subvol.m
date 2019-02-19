function fname = subvol(V,bb)
% Extract a subvolume
% FORMAT subvol(V,bb)
% V  - mapped file
% bb - bounding box
%
% Example:
%     V = spm_vol(spm_select(1,'image'));
%     subvol(V,[32 64 ; 1 64 ; 1 48]');
%
% %W% %E%

bb      = round(bb);
bb      = sort(bb);
bb(1,:) = max(bb(1,:),[1 1 1]);
bb(2,:) = min(bb(2,:),V.dim(1:3));

VO            = V;
[pth,nam,ext] = fileparts(V.fname);
fname         = fullfile(pth,['sv_' nam '.nii']);
VO.fname      = fname;
VO.dim(1:3)   = diff(bb)+1;
VO.mat        = V.mat*spm_matrix((bb(1,:)-1));

VO = spm_create_vol(VO);
for z=1:VO.dim(3),
    M   = V.mat\VO.mat*spm_matrix([0 0 z]);
    img = spm_slice_vol(V,M,VO.dim(1:2),0);
    VO  = spm_write_plane(VO,img,z);
end;

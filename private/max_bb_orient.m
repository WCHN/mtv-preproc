function [mat,dm] = max_bb_orient(Nii,vx)
% Calculate orientation matrix and dimensions from maximum bounding-box
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

mn = [ Inf  Inf  Inf]';
mx = [-Inf -Inf -Inf]';
for i=1:numel(Nii)
    N = numel(Nii{i});
    for n=1:N
        d = size(Nii{i}(n).dat);
        if numel(d)==2, d(3) = 1; end

        t = uint8(0:7);
        c = diag(d+1)*double([bitshift(bitand(t,bitshift(uint8(1),1-1)),1-1)
                              bitshift(bitand(t,bitshift(uint8(1),2-1)),1-2)
                              bitshift(bitand(t,bitshift(uint8(1),3-1)),1-3)]);
        c = bsxfun(@plus,Nii{i}(n).mat(1:3,1:3)*c,Nii{i}(n).mat(1:3,4));
        mx = max(mx,max(c,[],2));
        mn = min(mn,min(c,[],2));
    end
end
mat = spm_matrix(mn-1)*diag([vx 1])*spm_matrix(-[1 1 1]);
dm  = ceil((mat\[mx'+1 1]')');
dm  = dm(1:3);
%==========================================================================
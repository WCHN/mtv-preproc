function mat = rigidly_realign(mat,offset,rotation)
% Rigidly (translation and rotation) modify an orientation matrix
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

B                = zeros(4,4,6);
B(1,4,1)         = 1;
B(2,4,2)         = 1;
B(3,4,3)         = 1;
B([1,2],[1,2],4) = [0 1;-1 0];
B([3,1],[3,1],5) = [0 1;-1 0];
B([2,3],[2,3],6) = [0 1;-1 0];

Nr = size(B,3);
q  = zeros(1,Nr);

q(1:3) = offset(1:3);
q(4:6) = rotation(1:3);

R = spm_dexpm(q,B);

mat = R*mat;
%==========================================================================
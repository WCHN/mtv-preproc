function dat = init_dat(Nii,mat,dm,window,gap)
% Initialise projection matrices for super-resolution
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C   = numel(Nii);
dat = struct('mat',[],'dm',[],'N',[],'A',[]);
for c=1:C % Loop over channels
    if iscell(Nii(c))
        dat(c) = init_A(Nii{c},mat,dm,window{c},gap{c});         
    else
        dat(c) = init_A(Nii(c),mat,dm,window{c},gap{c});         
    end
end    
%==========================================================================

%==========================================================================
function dat = init_A(Nii,mat,dm,window,gap)
% Initialise projection matrices (stored in dat struct)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

N       = numel(Nii); % Number of LR images
dat.mat = mat;
dat.dm  = dm;   
dat.N   = N;
vs      = sqrt(sum(mat(1:3,1:3).^2));
for n=1:N % Loop over LR images
    mat_n = Nii(n).mat;    
    dm_n  = Nii(n).dat.dim;
    
    dat.A(n).mat = mat_n;    
    dat.A(n).dm  = dm_n;
    dat.A(n).win = window{n};
            
    M          = mat\mat_n;
%     R          = (M(1:3,1:3)/diag(sqrt(sum(M(1:3,1:3).^2))))';
%     dat.A(n).S = blur_fun(dm,R,sqrt(sum(M(1:3,1:3).^2)));
%     dat.A(n).S = blur_function(dm,M);

    S = sqrt(sum(M(1:3,1:3).^2));
    R = M(1:3,1:3)/diag(S);
    S = S - gap{n}./vs;
    M = R*diag(S);
    
    dat.A(n).J = single(reshape(M, [1 1 1 3 3]));
end
%==========================================================================

%==========================================================================
function f = blur_fun(dm,mat,vx)
if nargin<1, dm = [64 64]; end
if nargin<2, mat = eye(numel(dm)); end
if nargin<3, vx = ones(1,numel(dm)); end

if any(size(mat)~=numel(dm)) || numel(vx)~=numel(dm), error('Incompatible dimensions.'); end

% Grid in frequency space
r        = cell(1,numel(dm));
for i=1:numel(dm) 
    r{i} = single([0:ceil(dm(i)/2-1) -floor(dm(i)/2):-1]'*pi/dm(i)); 
end
X        = cell(1,numel(dm));
[X{:}]   = ndgrid(r{:});
clear r

% Transform
Y            = cell(size(X));
for i=1:numel(dm)
    Y{i}     = single(0);
    for j=1:numel(dm) 
        Y{i} = Y{i} + mat(i,j)*X{j}; 
    end
end
clear X

% Window function
f     = single(0);
for i=1:numel(dm) 
    f = f + Y{i}.^2; 
end    
f     = ((cos(min(f,pi^2/4)*4/pi) + 1)/2);

% Incorporate voxel size
for i=1:numel(dm)
    tmp                 = sin((vx(i))*Y{i})./(Y{i}.*cos(Y{i}/pi^(1/2)));
    tmp(~isfinite(tmp)) = vx(i);
    f                   = f.*tmp;
end
%==========================================================================
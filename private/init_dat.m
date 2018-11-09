function dat = init_dat(Nii,mat,dm)
% Initialise projection matrices for super-resolving MRIs
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C   = numel(Nii);
dat = struct('mat',[],'dm',[],'N',[],'A',[]);
for c=1:C        
    dat(c) = init_A(Nii(c),mat,dm);         
end    
%==========================================================================

%==========================================================================
function dat = init_A(Nii,mat,dm)
% Initialise projection matrices for super-resolving MRIs
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

N = numel(Nii); % Number of LR images of the same contrast  

dat.mat = mat;
dat.dm  = dm;   
dat.N   = N;
for n=1:N % Loop over LR images of the same contrast
    mat_n = Nii(n).mat;    
    dm_n  = Nii(n).dat.dim;
    
    dat.A(n).mat = mat_n;    
    dat.A(n).dm  = dm_n;
            
    M          = mat\mat_n;
    R          = (M(1:3,1:3)/diag(sqrt(sum(M(1:3,1:3).^2))))';
    dat.A(n).S = blur_fun(dm,R,sqrt(sum(M(1:3,1:3).^2)));
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
if 1
    f     = single(0);
    for i=1:numel(dm) 
        f = f + Y{i}.^2; 
    end    
    f     = ((cos(min(f,pi^2/4)*4/pi) + 1)/2);
else
    w_func     = @gausswin;
    sd         = 50;
    wr         = window(w_func,dm(2),sd);
    wc         = window(w_func,dm(1),sd);
    wz         = window(w_func,dm(3),sd);
    [wr,wc,wz] = meshgrid(wr,wc,wz);
    f          = wr.*wc.*wz;
    f          = single(abs(fftn(f)));    
end
% figure(111); imshow3D(f)

% Incorporate voxel size
for i=1:numel(dm)
    tmp                 = sin((vx(i))*Y{i})./(Y{i}.*cos(Y{i}/pi^(1/2)));
    tmp(~isfinite(tmp)) = vx(i);
    f                   = f.*tmp;
end
%==========================================================================
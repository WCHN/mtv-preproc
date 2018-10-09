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
function f = blur_fun(d,M,n)
if nargin<1, d = [64 64]; end
if nargin<2, M = eye(numel(d)); end
if nargin<3, n = ones(1,numel(d)); end

if any(size(M)~=numel(d)) || numel(n)~=numel(d), error('Incompatible dimensions.'); end

r    = cell(1,numel(d));
X    = cell(1,numel(d));
for i=1:numel(d), r{i} = single([0:ceil(d(i)/2-1) -floor(d(i)/2):-1]'*pi/d(i)); end
[X{:}] = ndgrid(r{:});

Y  = cell(size(X));
for i=1:numel(d)
    Y{i} = single(0);
    for j=1:numel(d), Y{i} = Y{i} + M(i,j)*X{j}; end
end
clear X

f  = single(0);
for i=1:numel(d), f = f + Y{i}.^2; end

f  = ((cos(min(f,pi^2/4)*4/pi)+1)/2);  % Some sort of window function

for i=1:numel(d)
    tmp = sin((n(i))*Y{i})./(Y{i}.*cos(Y{i}/pi^(1/2)));
    tmp(~isfinite(tmp)) = n(i);
    f = f.*tmp;
end
%==========================================================================
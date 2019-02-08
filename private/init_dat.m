function dat = init_dat(Nii,mat,dm,window,gap)
% Initialise projection matrices for super-resolution
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 4, window = 2; end
if nargin < 5, gap    = 0; end

% Get rigid basis
B = get_rigid_basis;

% Slice profile
window = get_window(window,Nii);

% Slice gap
gap = get_slice_gap(gap,Nii);

C   = numel(Nii);
dat = struct('mat',[],'dm',[],'N',[],'A',[]);
for c=1:C % Loop over channels
    if iscell(Nii(c))
        dat(c) = init_A(Nii{c},mat,dm,window{c},gap{c},B);         
    else
        dat(c) = init_A(Nii(c),mat,dm,window{c},gap{c},B);         
    end
end    
%==========================================================================

%==========================================================================
function dat = init_A(Nii,mat,dm,window,gap,B)
% Initialise projection matrices (stored in dat struct)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

N       = numel(Nii); % Number of LR images
Nq      = size(B,3);  % Number of affine parameters
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
    dat.A(n).gap = gap{n};
    dat.A(n).q   = zeros([Nq 1],'single'); 
    
    R = spm_dexpm(dat.A(n).q,B); % Rigid matrix
    M = mat\R*mat_n;
%     R          = (M(1:3,1:3)/diag(sqrt(sum(M(1:3,1:3).^2))))';
%     dat.A(n).S = blur_fun(dm,R,sqrt(sum(M(1:3,1:3).^2)));
%     dat.A(n).S = blur_function(dm,M);

    % Include slice-gap
    M = model_slice_gap(M,gap{n},vs);
    
    dat.A(n).J = single(reshape(M, [1 1 1 3 3]));
end
%==========================================================================

%==========================================================================
function gap = get_slice_gap(gap,Nii_x)
% Construct slice-gap
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(Nii_x);

if ~iscell(gap)
    gap = {gap};
end
gap = padarray(gap, [0 max(0,C-numel(gap))], 'replicate', 'post');
for c=1:C
    if ~iscell(gap{c})
        gap{c} = {gap{c}};
    end
    gap{c} = padarray(gap{c}, [0 max(0,numel(Nii_x{c})-numel(gap{c}))], 'replicate', 'post');
    for i=1:numel(gap{c})
        if isempty(gap{c}{i})
            gap{c}{i} = 2;
        end
        gap{c}{i} = padarray(gap{c}{i}, [0 max(0,3-numel(gap{c}{i}))], 'replicate', 'post');
    end
end
%==========================================================================

%==========================================================================
function window = get_window(window,Nii_x)
% Construct slice-profile
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(Nii_x);

if ~iscell(window)
    window = {window};
end
window = padarray(window, [0 max(0,C-numel(window))], 'replicate', 'post');
for c=1:C
    if ~iscell(window{c})
        window{c} = {window{c}};
    end
    window{c} = padarray(window{c}, [0 max(0,numel(Nii_x{c})-numel(window{c}))], 'replicate', 'post');
    for i=1:numel(window{c})
        if isempty(window{c}{i})
            window{c}{i} = 2;
        end
        window{c}{i} = padarray(window{c}{i}, [0 max(0,3-numel(window{c}{i}))], 'replicate', 'post');
    end
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
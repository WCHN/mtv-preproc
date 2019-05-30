function [sd,info] = noise_estimate_ct(Nii,K)
% Estimate noise in a CT image by fitting a variational Gaussian mixture
% model (GMM) to a range of its intensity histogram, and taking the
% component with a mean estimate closest to the value of air (HU = -1000)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 2, K = 2; end
   
% Get image date (vectorised)
f = Nii.dat(:);
f(~isfinite(f)) = [];
f(f == min(f))  = [];
f(f == max(f))  = [];

% Remove some voxels
mx = -990;
mn = -1023;
f(f >= mx)    = [];
f(f <  mn)    = [];

% Histogram bin voxels
x = min(f):max(f);
h = hist(f(:),x);

% Transpose
x = x';
h = h';

% Intensity distribution hyper-parameters
% MU0 = mean(f)*ones([1 K]) + K*randn([1 K]);
MU0 = [-1030 -1000];
b0  = ones([1 K]);
n0  = ones([1 K]);
V0  = ones([1 1 K]);

% Fit GMM
% ['kmeans'],'linspace','prior','sample','uniform',
% Start       = -1000*ones([1 K]) + K*randn([1 K]);
[~,MU,A,PI] = spm_gmm(x,K,h,'BinWidth',1,'GaussPrior',{MU0,b0,V0,n0}, ...
                      'Tolerance',1e-6,'Start','prior','Verbose',0,'Prune',true);

% Compute variance for each class
V = zeros(size(A));
for k=1:size(V,3)
    V(:,:,k) = inv(A(:,:,k));
end

% Get standard deviation of class closest to air (-1000)
[~,ix] = min(abs(abs(squeeze(MU)) - 1000));
sd     = sqrt(squeeze(V(:,:,ix)));

if nargout >= 2
    % Plot histogram + GMM fit
    p = zeros([numel(x) K],'single');
    for k=1:K
        p(:,k) = PI(k)*mvnpdf(x(:),MU(:,k),V(:,:,k));
    end
    sp = sum(p,2) + eps;    
    md = mean(diff(x));
    
    info    = struct;
    info.x  = x;
    info.f  = f;
    info.h  = h;
    info.p  = p;
    info.sp = sp;
    info.md = md;   
end
%==========================================================================
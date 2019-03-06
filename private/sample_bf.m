function bf = sample_bf(scl,lat,vs,fwhm,prm,Verbose)
% Sample a bias-field.
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 1, scl      = 1; end
if nargin < 2, lat      = [256 256 150]; end
if nargin < 3, vs       = [1 1 1]; end
if nargin < 4, fwhm     = 100; end
if nargin < 5, prm      = [1e-2 0 1e-3]; end
if nargin < 6, Verbose  = false; end

if scl == 0
    bf = ones(lat);
    return
end

nbcmp = fwhm2nbcomp(lat, vs, fwhm);

L     = regulariser(prm, lat, nbcmp, vs, 'n');
S     = inv(L);
[U,S] = svd(S);
S     = sqrt(S);
c     = reshape(U*S*randn(size(L,1),1), nbcmp);

[B1,B2,B3] = dcbasis(lat, nbcmp);

bf = zeros(lat);
for z=1:lat(3)
    bf(:,:,z) = exp(scl*transf(B1,B2,B3(z,:),c));
end

if Verbose || nargin == 0 || nargin == 1
    figure(111)
    imagesc3d(bf); colorbar
end
%==========================================================================

%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
return;
%==========================================================================

%==========================================================================
function nbcmp = fwhm2nbcomp(lattice, vs, fwhm)
% FORMAT nbcmp = spm_bias_lib('fwhm2nbcomp', lattice, voxel_size, fwhm)
%
% lattice      - Dimensions of the lattice [dx dy ...]
% voxel_size   - Voxel size of the lattice [vx vy ...]
% fwhm         - Full-width half-max of the highest frequency basis (mm)
%
% The number of components is chosen so that the full-width half-max of 
% the highest frequency basis function is smaller than fwhm. The bias 
% field cannot model effects whose spatial frequency is higher than this 
% value.
%
% If only one value is provided for voxel_size or fwhm, the same value is
% used along all dimensions.

% -------------------------------------------------------------------------
% Preprocess input arguments
ndim = numel(lattice);
vs = reshape(vs, 1, []);
if numel(vs) < ndim
    vs = padarray(vs, [0 ndim-numel(vs)], 'replicate', 'post');
end
fwhm = reshape(fwhm, 1, []);
if numel(fwhm) < ndim
    fwhm = padarray(fwhm, [0 ndim-numel(fwhm)], 'replicate', 'post');
end

% -------------------------------------------------------------------------
% Compute number of components per direction
nbcmp = ceil(2 * vs .* lattice ./ fwhm);
nbcmp = max(nbcmp, 1);
%==========================================================================

%==========================================================================
function varargout = dcbasis(lattice, nb_component)
% FORMAT [Bx,By,Bz,...] = spm_bias_lib('dcbasis', lattice, nb_component)
%
% lattice      - Dimensions of the lattice [dx dy ...]
% nb_component - Number of basis functions along each dimension [nx ny ...]
%
% Bx - Smooth basis along the x dimension [dx*nx] 
% By - Smooth basis along the y dimension [dy*ny]
% ...
%
% There are as many basis objects as elements in `lattice`

ndim = numel(lattice);

% -------------------------------------------------------------------------
% Preprocess input arguments
nb_component = reshape(nb_component, 1, []);
if numel(nb_component) < ndim
    nb_component = padarray(nb_component, [0 ndim-numel(nb_component)], 'replicate', 'post');
end

% -------------------------------------------------------------------------
% Compute each basis
varargout = cell(1,min(ndim, nargout));
for d=1:min(ndim, nargout)
    varargout{d} = spm_dctmtx(lattice(d),nb_component(d));
end
%==========================================================================

%==========================================================================
function L = regulariser(mode, lattice, nb_component, vs, bnd)
% FORMAT L = regulariser(param, lattice, nb_component, voxel_size)
% FORMAT L = regulariser(mode,  lattice, nb_component, voxel_size)
%
% param        - Parameters for absolute, membrane and bending energies
% mode         - Name of a single energy ('absolute'/'membrane'/'bending')
% lattice      - Dimensions of the lattice [dx dy ...]
% nb_component - Number of basis functions along each dimension [nx ny ...]
% voxel_size   - Voxel size of the lattice [vx vy ...]
%
% L            - Precision matrix [(nx*ny*...)^2]
%
% If numerical parameters are provided, a weighted combination of the  
% three types of regularisation is returned.
% If an energy name is provided, the matrix that allows to compute it is
% returned (without weighting: the regularisation parameter should be 
% multiplied with this matrix)
%
% If only one value is provided for nb_component or voxel_size, the
% same value is used along all dimensions.

if nargin < 5
    bnd = 'neumann';
end

% -------------------------------------------------------------------------
% Special case: mixture of regularisers
if ~ischar(mode)
    param = mode;
    L = 0;
    for i=1:numel(param)
        if param(i) ~= 0
            switch i
                case 1
                    mode = 'absolute';
                case 2
                    mode = 'membrane';
                case 3
                    mode = 'bending';
                case 4
                    mode = 'linearelastic1';
                case 5
                    mode = 'linearelastic2';
            end
            L1 = param(i) * regulariser(mode, lattice, nb_component, vs, bnd);
            if numel(L) == 1 || size(L,1) == size(L1,1)
                L = L + L1;
            else
                L0   = L;
                nprm = size(L,1);
                ndim = size(L1,1)/nprm;
                L = zeros(ndim*nprm);
                for d=1:ndim
                    L(((d-1)*nprm+1):d*nprm,((d-1)*nprm+1):d*nprm) = L0;
                end
                clear L0
                L = L + L1;
            end
        end
    end
    return
end


% -------------------------------------------------------------------------
% Preprocess input arguments
ndim = numel(lattice);
nb_component = reshape(nb_component, 1, []);
if numel(nb_component) < ndim
    nb_component = padarray(nb_component, [0 ndim-numel(nb_component)], 'replicate', 'post');
end
if nargin < 4
    vs = 1;
end
vs = reshape(vs, 1, []);
if numel(vs) < ndim
    vs = padarray(vs, [0 ndim-numel(vs)], 'replicate', 'post');
end

% -------------------------------------------------------------------------
% Mode-specific options
switch lower(mode)
    case {'absolute' 'abs' 'a'}
        maxdiff = 0;
    case {'membrane' 'mem' 'm' ...
          'linear-elastic1' 'linearelastic1' 'le1' ...
          'linear-elastic2' 'linearelastic2' 'le2'}
        maxdiff = 1;
    case {'bending' 'ben' 'b'}
        maxdiff = 2;
    otherwise
        error('Unknown mode %s, should be ''absolute'', ''membrane'' or ''bending''.', mode);
end

% -------------------------------------------------------------------------
% Compute each basis + square it
switch lower(bnd)
    case {0, 'circulant', 'circ', 'c'}
        mtxfun = @spm_dftmtx;
    case {1, 'neumann', 'neu', 'n'}
        mtxfun = @spm_dctmtx;
    case {2, 'dirichlet', 'dir', 'd'}
        mtxfun = @spm_dstmtx;
    otherwise
        error('Unknown boundary condition');
end

basis = cell(ndim, maxdiff + 1);
nbprm = 1;
for d=1:ndim
    for diff=0:maxdiff
        switch diff
            case 0
                basis{d,diff+1} = mtxfun(lattice(d),nb_component(d));
                nbprm = nbprm * size(basis{d,diff+1}, 2);
            case 1
                basis{d,diff+1} = mtxfun(lattice(d),nb_component(d),'diff') / vs(d);
            case 2
                basis{d,diff+1} = mtxfun(lattice(d),nb_component(d),'diff2') / vs(d)^2;
        end
        if any(strcmpi(mode, {'absolute' 'abs' 'a' 'membrane' 'mem' 'm' 'bending' 'ben' 'b'}))
            basis{d,diff+1} = basis{d,diff+1}.' * basis{d,diff+1};
        end
    end
end

% -------------------------------------------------------------------------
% Compute precision matrix
switch lower(mode)
    case {'absolute' 'abs' 'a'}
        L = 1;
        for d=1:ndim
            L = spm_krutil(basis{d,1}, L);
        end
    case {'membrane' 'mem' 'm'}
        L = 0;
        for dd=1:ndim               % Which dimension to differentiate
            L1 = 1;
            for d=1:ndim            % Kronecker loop
                if d == dd
                    L1 = spm_krutil(basis{d,2}, L1);
                else
                    L1 = spm_krutil(basis{d,1}, L1);
                end
            end
            L = L + L1;
        end
    case {'bending' 'ben' 'b'}
        L = 0;
        for dd1=1:ndim              % First dimension to differentiate
            L1 = 1;
            for d=1:ndim            % Kronecker loop
                if d == dd1
                    L1 = spm_krutil(basis{d,3}, L1);
                else
                    L1 = spm_krutil(basis{d,1}, L1);
                end
            end
            L = L + L1;
            for dd2=dd1+1:ndim      % Second dimension to differentiate
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == dd1 || d == dd2
                        L1 = spm_krutil(basis{d,2}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}, L1);
                    end
                end
                L = L + 2 * L1;
            end
        end
    case {'linear-elastic1' 'linearelastic1' 'le1'}
        L = zeros(nbprm,ndim,nbprm,ndim);
        for h1=1:ndim               % First Hessian dimension
            for dd=1:ndim          % First dimension to differentiate
                if dd == h1
                    coeff = 1;
                else
                    coeff = 0.5;
                end
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == dd
                        L1 = spm_krutil(basis{d,2}.' * basis{d,2}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                    end
                end
                L(:,h1,:,h1) = L(:,h1,:,h1) + coeff * reshape(L1, [nbprm 1 nbprm]);
            end
            
            for h2=h1+1:ndim        % Second Hessian dimension
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == h1
                        L1 = spm_krutil(basis{d,1}.' * basis{d,2}, L1);
                    elseif d == h2
                        L1 = spm_krutil(basis{d,2}.' * basis{d,1}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                    end
                end
                L(:,h1,:,h2) = L(:,h1,:,h2) + 0.5 * reshape(L1,  [nbprm 1 nbprm]);
                L(:,h2,:,h1) = L(:,h2,:,h1) + 0.5 * reshape(L1', [nbprm 1 nbprm]);
            end
            
        end
        L = reshape(L, nbprm*ndim, nbprm*ndim);
    case {'linear-elastic2' 'linearelastic2' 'le2'}
        L = zeros(nbprm,ndim,nbprm,ndim);
        for h1=1:ndim               % First Hessian dimension
            L1 = 1;
            for d=1:ndim        % Kronecker loop
                if d == h1
                    L1 = spm_krutil(basis{d,2}.' * basis{d,2}, L1);
                else
                    L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                end
            end
            L(:,h1,:,h1) = L(:,h1,:,h1) + 0.5 * reshape(L1, [nbprm 1 nbprm]);
            
            for h2=h1+1:ndim        % Second Hessian dimension
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == h1
                        L1 = spm_krutil(basis{d,2}.' * basis{d,1}, L1);
                    elseif d == h2
                        L1 = spm_krutil(basis{d,1}.' * basis{d,2}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                    end
                end
                L(:,h1,:,h2) = L(:,h1,:,h2) + 0.5 * reshape(L1,  [nbprm 1 nbprm]);
                L(:,h2,:,h1) = L(:,h2,:,h1) + 0.5 * reshape(L1', [nbprm 1 nbprm]);
            end
            
        end
        L = reshape(L, nbprm*ndim, nbprm*ndim);
end
%==========================================================================
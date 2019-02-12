function [dat,ll1,armijo] = update_rigid(Nii_x,Nii_y,dat,tau,armijo,num_workers,p)
% Optimise the rigid alignment between observed images and their 
% corresponding channel's reconstruction. This routine runs in parallel
% over image channels. The observed lowres images are here denoted by f and
% the reconstructions by mu.
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Some parameters from options struct
meancrct = p.Results.MeanCorrectRigid;
nitgn    = p.Results.IterGaussNewtonRigid; 
speak    = p.Results.Verbose; 

% More parameters
B  = get_rigid_basis; % Get rigid basis
C  = numel(dat);      % Number of channels
Nq = size(B,3);       % Number of registration parameters

if num_workers > 0
    speak = min(speak,1);
end

%--------------------------------------------------------------------------
% Start updating, for each observation
%--------------------------------------------------------------------------

ll1 = zeros(1,C);
% for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers) % Loop over channels
    
    % Boundary used to model HR image  
    spm_field('boundary',1);
    pushpull('boundary',1); 
    
    [dat(c),ll1(c),armijo{c}] = update_channel(Nii_x(c),Nii_y(c),dat(c),B,tau{c},armijo{c},speak,nitgn,c);    
    
end % End loop over channels

if meancrct
    %----------------------------------------------------------------------
    % Mean correct the rigid-body transforms
    %----------------------------------------------------------------------

    % Sum all rigid body parameters
    sq  = zeros([Nq 1]);
    cnt = 0;
    for c=1:C
        N = numel(dat(c).A);
        for n=1:N
            sq  = sq + dat(c).A(n).q;
            cnt = cnt + 1;
        end
    end

    % Compute mean
    q_avg = sq/cnt;

    % Update all q and J
    for c=1:C

        Mmu  = dat(c).mat;        
        vsmu = sqrt(sum(Mmu(1:3,1:3).^2));

        N = numel(dat(c).A);
        for n=1:N

            Mf = dat(c).A(n).mat;                   
            q  = dat(c).A(n).q;

            q = q - q_avg; 

            R  = spm_dexpm(q,B);
            M  = Mmu\R*Mf;
            Mg = model_slice_gap(M,dat(c).A(n).gap,vsmu);    
            J  = single(reshape(Mg, [1 1 1 3 3]));

            % Update dat struct
            dat(c).A(n).q = q;        
            dat(c).A(n).J = J;
        end  
    end
end
%==========================================================================

%==========================================================================
function [dat,sll,armijo] = update_channel(Nii_x,Nii_y,dat,B,tau,armijo,speak,nitgn,c)

% Parameters
Nq          = size(B,3);             % Number of registration parameters
lkp         = [1 4 5; 4 2 6; 5 6 3]; % Que?
nlinesearch = 4;                     % Number of line-searches    
method      = dat.method;

% Cell-array of observations    
f = get_nii(Nii_x); 
N = numel(f); % Number of observations

% Template parameters
mu   = get_nii(Nii_y);
Mmu  = dat.mat;        
vsmu = sqrt(sum(Mmu(1:3,1:3).^2));

sll = 0;
for n=1:N % Loop over observed images (of channel c)

    % Observation parameters
    dmf  = size(f{n});        
    Mf   = dat.A(n).mat;            
    tauf = tau(n);        % Noise precision

    for gnit=1:nitgn % Loop over Gauss-Newton iterations
        
        oq = dat.A(n).q;                        
        oJ = dat.A(n).J;
        
        %------------------------------------------------------------------
        % Compute gradient and Hessian (slice-wise)
        %------------------------------------------------------------------

        % Differentiate Rq w.r.t. q, store in dRq
        [R,dR]         = spm_dexpm(oq,B);                     
        dRq            = zeros(4,4,Nq);
        for i=1:Nq
            dRq(:,:,i) = Mmu\dR(:,:,i)*Mf;
        end

        % Make interpolation grids: yf - image, ymu - template-2-image
        M     = Mmu\R*Mf;
        yf    = get_id(dmf);        
        ymu2f = affine_transf(M,yf);

        % Build gradient and Hessian
        g  = zeros([Nq 1], 'single');
        H  = zeros([Nq Nq],'single');
        ll = 0;
        for z=1:dmf(3) % Loop over slices

            % Compute matching-term part (log likelihood)
            [llz,gz,Hz] = meansq_objfun_slice(f{n},mu,ymu2f,dat.A(n),tauf,speak,z,method);
            ll          = ll + llz;           

            % Add dRq to gradients
            dAff = cell(Nq,3);
            for i=1:Nq
                for d=1:3
                    tmp       = dRq(d,1,i)*yf(:,:,z,1) + dRq(d,2,i)*yf(:,:,z,2) + dRq(d,3,i)*yf(:,:,z,3) + dRq(d,4,i);
                    dAff{i,d} = tmp(:);
                end
            end
            for d=1:3
                tmp = gz(:,:,d);
                tmp = tmp(:)';
                for i=1:Nq
                    g(i) = g(i) + tmp*dAff{i,d};
                end
            end

            % Add dRq to Hessian
            for d1=1:3
                for d2=1:3
                    tmp1 = Hz(:,:,lkp(d1,d2));
                    tmp1 = tmp1(:);
                    for i1=1:Nq
                        tmp2 = (tmp1.*dAff{i1,d1})';
                        for i2=i1:Nq % Only compute upper/lower triangle of symmetric matrix
                            H(i1,i2) = H(i1,i2) + tmp2*dAff{i2,d2};
                        end
                    end
                end
            end

        end % End loop over slices

        if speak >= 1
            fprintf('   | c=%i, n=%i, gn=%i, ls=0 | ll=%g\n', c, n, gnit, ll); 
        end

        % Fill in missing triangle
        for i1=1:Nq
            for i2=1:Nq
                H(i2,i1) = H(i1,i2);
            end
        end                

        %------------------------------------------------------------------
        % Update q by Gauss-Newton optimisation
        %------------------------------------------------------------------

        % Compute update step from gradient and Hessian
        H      = loaddiag(H);
        Update = H\g;

        % Start line-search                       
        oll = ll;        
        for linesearch=1:nlinesearch

            % Take step
            q = oq - armijo(n)*Update;

            % Compute new parameters
            R  = spm_dexpm(q,B);
            M  = Mmu\R*Mf;
            Mg = model_slice_gap(M,dat.A(n).gap,vsmu);    
            J  = single(reshape(Mg, [1 1 1 3 3]));

            % Update dat struct
            dat.A(n).q = q;        
            dat.A(n).J = J;

            % Compute new log-likelihood
            ymu2f = affine_transf(M,yf);
            ll    = meansq_objfun(f{n},mu,ymu2f,dat.A(n),tauf,speak,method,true);

            if ll > oll
                if speak >= 1
                    fprintf('   | c=%i, n=%i, gn=%i, ls=%i | ll=%g | :o) | q=%s\n', c, n, gnit, linesearch, ll, sprintf(' %2.3f', q)); 
                end

                armijo(n) = min(1.2*armijo(n),1);

                break;
            else
                if speak >= 1
                    fprintf('   | c=%i, n=%i, gn=%i, ls=%i | ll=%g | :''(\n', c, n, gnit, linesearch, ll); 
                end
                
                % Revert to previous values in dat struct
                dat.A(n).q = oq;        
                dat.A(n).J = oJ;
                
                if armijo(n) < eps('single')
                    % We are probably close to the optimum, so cancel
                    % line-search
                    break
                end
                
                armijo(n) = 0.5*armijo(n);
            end
        end

        if ll <= oll
            % Use old neg log-likelihood
            ll = oll;
        end
        
    end % End loop over Gauss-Newton iterations

    sll = sll + ll;
    
end % End loop over observations
%==========================================================================
    
%==========================================================================
function ll = meansq_objfun(f,mu,y,An,tau,speak,method,show_moved)
if nargin < 8, show_moved = 0; end

dm = size(f);
ll = 0;
for z=1:dm(3)
    ll = ll + meansq_objfun_slice(f,mu,y,An,tau,speak,z,method,show_moved);
end
%==========================================================================

%==========================================================================
function [ll,g,H] = meansq_objfun_slice(f,mu,y,An,tau,speak,z,method,show_moved)
if nargin < 9, show_moved = 0; end

dm                = size(f); % Observation dimensions
mu(~isfinite(mu)) = 0; 

% Move template to image space
dmu = cell(1,3);
if nargout >= 2
    if strcmp(method,'superres')
        [mu,dmu{1},dmu{2},dmu{3}] = pushpull('pull',single(mu),single(y(:,:,z,:)),single(An.J),double(An.win)); 
    elseif strcmp(method,'denoise')
        [~,dmu{1}]     = spm_diffeo('bsplins',mu,y(:,:,z,:),[2 0 0  0 0 0]);        
        [~,~,dmu{2}]   = spm_diffeo('bsplins',mu,y(:,:,z,:),[0 2 0  0 0 0]);
        [~,~,~,dmu{3}] = spm_diffeo('bsplins',mu,y(:,:,z,:),[0 0 2  0 0 0]);
        
        mu = spm_diffeo('pull',mu,y(:,:,z,:));        
    end
else
    if strcmp(method,'superres')
        mu = pushpull('pull',single(mu),single(y(:,:,z,:)),single(An.J),double(An.win));    
    elseif strcmp(method,'denoise')
        mu = spm_diffeo('pull',mu,y(:,:,z,:));
    end
end
mu(~isfinite(mu)) = 0; 

if speak >= 2 && z == round((dm(3) + 1)/2)     
    % Some verbose    
    show_reg(mu,f(:,:,z),dmu,show_moved);
end

if nargout == 0, return; end

% Compute log-likelihood of slice
msk  = isfinite(f(:,:,z)) & f(:,:,z) ~= 0;
msk  = msk(:);
ftmp = f(:,:,z);
ll   = -0.5*tau*sum((double(ftmp(msk)) - double(mu(msk))).^2);

if nargout >= 2
    % Compute gradient    
    g = zeros([dm(1:2),3],'single');
    
    diff1        = mu - f(:,:,z);
    for d=1:3
        g(:,:,d) = diff1.*dmu{d};
    end
    
    if nargout >= 3
        % Compute Hessian
        H = zeros([dm(1:2),6],'single');
        
        H(:,:,1) = dmu{1}.*dmu{1};
        H(:,:,2) = dmu{2}.*dmu{2};
        H(:,:,3) = dmu{3}.*dmu{3};
        H(:,:,4) = dmu{1}.*dmu{2};
        H(:,:,5) = dmu{1}.*dmu{3};
        H(:,:,6) = dmu{2}.*dmu{3};                
    end
end

% Remove missing data from derivatives (represented by NaNs)
if nargout >= 2
    g(~isfinite(g)) = 0;
    if nargout >= 3
        H(~isfinite(H)) = 0;
    end
end
%==========================================================================

%==========================================================================
function y1 = affine_transf(Affine,y)
y1 = cat(4, Affine(1,1)*y(:,:,:,1) + Affine(1,2)*y(:,:,:,2) + Affine(1,3)*y(:,:,:,3) + Affine(1,4),...
            Affine(2,1)*y(:,:,:,1) + Affine(2,2)*y(:,:,:,2) + Affine(2,3)*y(:,:,:,3) + Affine(2,4),...
            Affine(3,1)*y(:,:,:,1) + Affine(3,2)*y(:,:,:,2) + Affine(3,3)*y(:,:,:,3) + Affine(3,4));
%==========================================================================            

%==========================================================================
function id = get_id(dm)
id        = cell(3,1);
[id{1:3}] = ndgrid(single(1:dm(1)),single(1:dm(2)),single(1:dm(3)));
id        = cat(4,id{:});
%==========================================================================

%==========================================================================
function A = loaddiag(A)
factor = 1e-7;
while rcond(A) < 1e-5
    A = A + factor * max([diag(A); eps]) * eye(size(A));
    factor = 10 * factor;
end
%==========================================================================

%==========================================================================
function show_reg(mu,f,dmu,show_moved)
if nargin < 4, show_moved = 0; end

figname = '(SPM) Rigid registration';
fig     = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);  

mxf = max(f(:));
if ~show_moved
%     clf(fig)
    subplot(2,3,1); imagesc(f', [0 mxf]); axis xy off; title('f');
    subplot(2,3,2); imagesc(mu',[0 mxf]); axis xy off; title('mu');
    subplot(2,3,3); imagesc(mu',[0 mxf]); axis xy off; title('nmu');
    subplot(2,3,4); imagesc(dmu{1}'); axis xy off; title('dmux');
    subplot(2,3,5); imagesc(dmu{2}'); axis xy off; title('dmuy');
    subplot(2,3,6); imagesc(dmu{3}'); axis xy off; title('dmuz');
else
    subplot(2,3,3); imagesc(mu',[0 mxf]); axis xy off; title('nmu');
end
drawnow
%==========================================================================    
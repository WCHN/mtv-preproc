function [dat,ll1] = update_rigid(Nii,dat,tau,num_workers,p)
% Update q with a Gauss-Newton algorithm
%
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
is3d = dat(1).dm(3) > 1;
B    = get_rigid_basis(is3d); % Get rigid basis
C    = numel(dat);            % Number of channels
Nq   = size(B,3);             % Number of registration parameters

if num_workers > 0
    speak = min(speak,1);
end

%--------------------------------------------------------------------------
% Start updating, for each observation
%--------------------------------------------------------------------------

ll1 = zeros(1,C);
% for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers) % Loop over channels
    
    set_boundary_conditions;
    
    [dat(c),ll1(c)] = update_channel(Nii.x(c),Nii.y(c),dat(c),B,tau{c},speak,nitgn,c);    
    
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
function [dat,sll] = update_channel(Nii_x,Nii_y,dat,B,tau,speak,nitgn,c)

% Parameters
Nq          = size(B,3);             % Number of registration parameters
lkp         = [1 4 5; 4 2 6; 5 6 3]; % Que?
nlinesearch = 6;                     % Number of line-searches    
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
    dmf  = [dmf 1];
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

        % Make interpolation grids
        M = Mmu\R*Mf;
        x = get_id(dmf);        
        y = affine_transf(M,x);

        % Build gradient and Hessian
        g  = zeros([Nq 1]);
        H  = zeros([Nq Nq]);
        ll = 0;
        for z=1:dmf(3) % Loop over slices

            % Compute matching-term part (log likelihood)
            [llz,gz,Hz] = meansq_objfun_slice(f{n},mu,y,dat.A(n),tauf,speak,z,method);
            ll          = ll + llz;           

            % Add dRq to gradients
            dAff = cell(Nq,3);
            for i=1:Nq
                for d=1:3
                    tmp       = dRq(d,1,i)*double(x(:,:,z,1)) + ...
                                dRq(d,2,i)*double(x(:,:,z,2)) + ...
                                dRq(d,3,i)*double(x(:,:,z,3)) + ...
                                dRq(d,4,i);
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

        if 0
            compare_numerical_derivative(g,f{n},mu,dat.A(n),x,tauf,Mmu,Mf,speak,method);
        end
        
        %------------------------------------------------------------------
        % Update q by Gauss-Newton optimisation
        %------------------------------------------------------------------

        % Compute update step from gradient and Hessian
        Update = H\g;

        % Start line-search                       
        oll    = ll;        
        armijo = 1;
        for linesearch=1:nlinesearch

            % Take step
            q = oq - armijo*Update;

            % Compute new parameters
            R  = spm_dexpm(q,B);
            M  = Mmu\R*Mf;
            Mg = model_slice_gap(M,dat.A(n).gap,vsmu);    
            J  = single(reshape(Mg, [1 1 1 3 3]));

            % Update dat struct
            dat.A(n).q = q;        
            dat.A(n).J = J;

            % Compute new log-likelihood
            y  = affine_transf(M,x);
            ll = meansq_objfun(f{n},mu,y,dat.A(n),tauf,speak,method,true);

            if ll > oll
                if speak >= 1
                    fprintf('   | c=%i, n=%i, gn=%i, ls=%i | ll=%g | :o) | q=%s\n', c, n, gnit, linesearch, ll, sprintf(' %2.3f', q)); 
                end
                
                break;
            else                
                % Revert to previous values in dat struct
                dat.A(n).q = oq;        
                dat.A(n).J = oJ;
                
                armijo = 0.5*armijo;
            end
        end

        if ll <= oll
            % Use old log-likelihood
            ll = oll;
            
            if speak >= 1
                fprintf('   | c=%i, n=%i, gn=%i, ls=%i | ll=%g | :''(\n', c, n, gnit, linesearch, ll); 
            end
        end
        
    end % End loop over Gauss-Newton iterations

    sll = sll + ll;
    
end % End loop over observations
%==========================================================================
    
%==========================================================================
function ll = meansq_objfun(f,mu,y,A,tau,speak,method,show_moved)
if nargin < 8, show_moved = 0; end

dm = size(f);
dm = [dm 1];
ll = 0;
for z=1:dm(3)
    ll = ll + meansq_objfun_slice(f,mu,y,A,tau,speak,z,method,show_moved);
end
%==========================================================================

%==========================================================================
function [ll,g,H] = meansq_objfun_slice(f,mu,y,A,tau,speak,z,method,show_moved)
if nargin < 9, show_moved = 0; end

dm                = size(f); % Observation dimensions
dm                = [dm 1];
mu(~isfinite(mu)) = 0; 

% Move template to image space
dmu = cell(1,3);
if nargout >= 2
    if strcmp(method,'superres')
        [mu,dmu{1},dmu{2},dmu{3}] = pushpull('pull',single(mu),single(y(:,:,z,:)),single(A.J),double(A.win)); 
    elseif strcmp(method,'denoise')               
        [mu,dmu{1},dmu{2},dmu{3}] = spm_diffeo('bsplins',mu,y(:,:,z,:),[1 1 1  0 0 0]); 
    end
else
    if strcmp(method,'superres')
        mu = pushpull('pull',single(mu),single(y(:,:,z,:)),single(A.J),double(A.win));    
    elseif strcmp(method,'denoise')
        mu = spm_diffeo('bsplins',mu,y(:,:,z,:),[1 1 1  0 0 0]); 
    end
end
mu(~isfinite(mu)) = 0; 

if speak >= 2 && z == floor((dm(3) + 1)/2)     
    % Some verbose    
    show_reg(mu,f(:,:,z),dmu,show_moved);
end

if nargout == 0, return; end

% Double
mu     = double(mu);
dmu{1} = double(dmu{1});
dmu{2} = double(dmu{2});
dmu{3} = double(dmu{3});

% Compute log-likelihood of slice
msk  = isfinite(f(:,:,z)) & f(:,:,z) ~= 0;
msk  = msk(:);
ftmp = f(:,:,z);
ll   = -0.5*tau*sum((double(ftmp(msk)) - mu(msk)).^2);

if nargout >= 2
    % Compute gradient    
    g = zeros([dm(1:2),3]);
    
    diff1        = mu - double(f(:,:,z));
    for d=1:3
        g(:,:,d) = diff1.*dmu{d};
    end
    
    if nargout >= 3
        % Compute Hessian
        H = zeros([dm(1:2),6]);
        
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
function y = affine_transf(Affine,x)
dm = size(x);
y  = cat(4, Affine(1,1)*x(:,:,:,1) + Affine(1,2)*x(:,:,:,2) + Affine(1,3)*x(:,:,:,3) + Affine(1,4),...
            Affine(2,1)*x(:,:,:,1) + Affine(2,2)*x(:,:,:,2) + Affine(2,3)*x(:,:,:,3) + Affine(2,4),...
            Affine(3,1)*x(:,:,:,1) + Affine(3,2)*x(:,:,:,2) + Affine(3,3)*x(:,:,:,3) + Affine(3,4));
if dm(3) == 1
    y(:,:,:,end) = 1;
end        
%==========================================================================            

%==========================================================================
function x = get_id(dm)
x        = cell(3,1);
[x{1:3}] = ndgrid(single(1:dm(1)),single(1:dm(2)),single(1:dm(3)));
x        = cat(4,x{:});
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

%==========================================================================
function compare_numerical_derivative(g,f,mu,A,x,tau,Mmu,Mf,speak,method)

% Parameters
dm   = size(x);
B    = get_rigid_basis(dm(3) > 1);               
oq   = A.q;   
Nq   = numel(oq);
vsmu = sqrt(sum(Mmu(1:3,1:3).^2));

% Modulate gradient w tau
g = tau*g;

% Get deformation
% R = spm_dexpm(oq,B);  
% M = Mmu\R*Mf;       
% y = affine_transf(M,x);

% Compute f(x)
% fx = meansq_objfun(f,mu,y,A,tau,speak,method);
% fx = -fx; % Because we differentiate the neg ll

%--------------------------------------------------------------------------
% Check derivatives numerically, by computing (f(x + h) - f(x))/h for
% various values of h, and by perturbing all individual parameters
%--------------------------------------------------------------------------

h = 10.^(-7:1:-1); % step-size

for n=1:Nq % Loop over rigid parameters
    
    fprintf('---------------------------\n')
    fprintf(' g(q(%i)) = %7.7f\n',n,g(n));
    
    for i=1:numel(h) % Loop over step sizes

        q    = oq;
        q(n) = q(n) + h(i);
        A.q  = q;
                               
        % Get deformation
        R = spm_dexpm(q,B);  
        M = Mmu\R*Mf;       
        y = affine_transf(M,x);
        
        Mg  = model_slice_gap(M,A.gap,vsmu);    
        J   = single(reshape(Mg, [1 1 1 3 3]));
        A.J = J;
        
        fxph = meansq_objfun(f,mu,y,A,tau,speak,method);
        fxph = -fxph; % Because we differentiate the neg ll

        q    = oq;
        q(n) = q(n) - h(i);
        A.q  = q;
                               
        % Get deformation
        R = spm_dexpm(q,B);  
        M = Mmu\R*Mf;       
        y = affine_transf(M,x);
        
        Mg  = model_slice_gap(M,A.gap,vsmu);    
        J   = single(reshape(Mg, [1 1 1 3 3]));
        A.J = J;
        
        fxmh = meansq_objfun(f,mu,y,A,tau,speak,method);
        fxmh = -fxmh; % Because we differentiate the neg ll
        
        gn = (fxph - fxmh)/(2*h(i));
        
        fprintf('ng(q(%i)) = %7.7f, for h = %g\n',n,gn,h(i));
    end
end
%==========================================================================
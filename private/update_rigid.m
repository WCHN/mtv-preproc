function [dat,armijo,ll1] = update_rigid(Nii,dat,tau,armijo,num_workers,p)
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
    
    [dat(c),ll1(c),armijo{c}] = update_channel(Nii.x(c),Nii.y(c),dat(c),B,tau{c},speak,nitgn,c,armijo{c});    
    
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
            dat(c).A(n).R = R;
        end  
    end
end
%==========================================================================

%==========================================================================
function [dat,sll,armijo] = update_channel(Nii_x,Nii_y,dat,B,tau,speak,nitgn,c,armijo)

% Parameters
Nq          = size(B,3);             % Number of registration parameters
lkp         = [1 4 5; 4 2 6; 5 6 3]; % Que?
nlinesearch = 12;                    % Number of line-searches    
method      = dat.method;

% Cell-array of observations    
f = get_nii(Nii_x); 
N = numel(f); % Number of observations

% Template parameters
mu    = get_nii(Nii_y);
Mmu   = dat.mat; 
dmmu  = [size(mu) 1];
dmmu  = dmmu(1:3);
vsmu  = sqrt(sum(Mmu(1:3,1:3).^2));
fovmu = abs(Mmu*[dmmu'+0.5; 1]-Mmu*[0.5*ones(3,1); 1]);
fovmu = fovmu(1:3);
is3d  = dmmu(3) > 1;

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
        oR = dat.A(n).R;
        
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
        g   = zeros([Nq 1]);
        H   = zeros([Nq Nq]);
        ll  = 0;
        onm = 0;
        for z=1:dmf(3) % Loop over slices

            % Compute matching-term part (log likelihood)
            [llz,dnm,gz,Hz] = meansq_objfun_slice(f{n},mu,y,dat.A(n),tauf,speak,z,method);
            ll              = ll + llz;           
            onm             = onm + dnm;
            
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

        % Fill in missing triangle
        for i1=1:Nq
            for i2=(i1+1):Nq
                H(i2,i1) = H(i1,i2);
            end
        end                

        % Regularise diagonal of Hessian a little bit        
        H = H + 1e-5*max(diag(H))*eye(size(H));
        
        if 0
            compare_numerical_derivatives(g,H,f{n},mu,dat.A(n),x,tauf,Mmu,Mf,speak,method);
        end
        
        %------------------------------------------------------------------
        % Update q by Gauss-Newton optimisation
        %------------------------------------------------------------------

        % Compute update step from gradient and Hessian
        Update = H\g;

        % Start line-search    
        armijo(n) = max_armijo(oq,Update,fovmu,is3d);
        oll       = ll;    
        onm       = nm;
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
            dat.A(n).R = R;

            % Compute new log-likelihood
            y       = affine_transf(M,x);
            [ll,nm] = meansq_objfun(f{n},mu,y,dat.A(n),tauf,speak,method,true);
            
            if ll/nm > oll/onm % && q_constraint(q,is3d)
                % Log-likelihood improved
                armijo(n) = min(1.2*armijo(n),1);

                break;
            else                
                % Revert to previous values in dat struct
                dat.A(n).q = oq;        
                dat.A(n).J = oJ;
                dat.A(n).R = oR;
                
                armijo(n) = 0.5*armijo(n);
            end
        end % End loop over line-search

        if ll/nm <= oll/onm
            % Use old log-likelihood
            ll = oll;                        
            nm = onm;
            
            if speak >= 1
                fprintf('   | c=%i, n=%i, gn=%i, ls=%i | ll=%10.1f, nm=%d | a=%7.5f, q=%s | :''(\n', c, n, gnit, linesearch, ll, nm0, armijo(n), sprintf(' %5.2f', dat.A(n).q)); 
            end
            
            % We didn't improve the objective function, so exit GN loop
            break
        elseif speak >= 1
            fprintf('   | c=%i, n=%i, gn=%i, ls=%i | ll=%10.1f, nm=%d  | a=%7.5f, q=%s | :o)\n', c, n, gnit, linesearch, ll, nm0, armijo(n), sprintf(' %5.2f', dat.A(n).q)); 
        end
        
    end % End loop over Gauss-Newton iterations

    sll = sll + ll;
    
end % End loop over observations
%==========================================================================
    
%==========================================================================
function [ll,nm] = meansq_objfun(f,mu,y,A,tau,speak,method,show_moved)
if nargin < 8, show_moved = 0; end

dm = size(f);
dm = [dm 1];
ll = 0;
nm = 0;
for z=1:dm(3)
    [dll,dnm] = meansq_objfun_slice(f,mu,y,A,tau,speak,z,method,show_moved);
    ll        = ll + dll;
    nm        = nm + dnm;
end
%==========================================================================

%==========================================================================
function [ll,nm,g,H] = meansq_objfun_slice(f,mu,y,A,tau,speak,z,method,show_moved)
if nargin < 9, show_moved = 0; end

dm = size(f); % Observation dimensions
dm = [dm 1];

% Move template to image space
dmu = cell(1,3);
if nargout >= 3
    if strcmp(method,'superres')
        [mu,dmu{1},dmu{2},dmu{3}] = pushpull('pull',mu,y(:,:,z,:),single(A.J),double(A.win)); 
    elseif strcmp(method,'denoise')
        % (JA) Not convinced by this, and would prefer the commented out version:
        [~,dmu{1},~,~] = spm_diffeo('bsplins',mu,y(:,:,z,:),[2 0 0  0 0 0]);
        [~,~,dmu{2},~] = spm_diffeo('bsplins',mu,y(:,:,z,:),[0 2 0  0 0 0]);
        [~,~,~,dmu{3}] = spm_diffeo('bsplins',mu,y(:,:,z,:),[0 0 2  0 0 0]);
        mu             = spm_diffeo('pull',mu,y(:,:,z,:));
%       [mu,dmu{1},dmu{2},dmu{3}] = spm_diffeo('bsplins',mu,y(:,:,z,:),[1 1 1  0 0 0]); 
    end
else
    if strcmp(method,'superres')
        mu = pushpull('pull',mu,y(:,:,z,:),single(A.J),double(A.win));    
    elseif strcmp(method,'denoise')
        mu = spm_diffeo('pull',mu,y(:,:,z,:));
%         mu = spm_diffeo('bsplins',mu,y(:,:,z,:),[1 1 1  0 0 0]);         
    end
end

% nm = sum(~isfinite(mu(:)));
% if nm > 0
%     disp(num2str(nm))
% end

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
msk  = get_msk(f(:,:,z),mu);
ftmp = f(:,:,z);
% ll   = -0.5*tau*sum((double(ftmp(msk)) - mu(msk)).^2);
ll   = -0.5*tau*sum(double(mu(msk).^2 - 2*mu(msk).*ftmp(msk)));

nm  = sum(msk);
msk = reshape(msk, dm(1:2));

if nargout >= 3
    % Compute gradient    
    g = zeros([dm(1:2),3]);
    
    diff1        = mu - double(f(:,:,z));
    for d=1:3
        g(:,:,d) = diff1.*dmu{d}.*msk;
    end
    
    if nargout >= 4
        % Compute Hessian
        H = zeros([dm(1:2),6]);
        
        H(:,:,1) = dmu{1}.*dmu{1}.*msk;
        H(:,:,2) = dmu{2}.*dmu{2}.*msk;
        H(:,:,3) = dmu{3}.*dmu{3}.*msk;
        H(:,:,4) = dmu{1}.*dmu{2}.*msk;
        H(:,:,5) = dmu{1}.*dmu{3}.*msk;
        H(:,:,6) = dmu{2}.*dmu{3}.*msk;             
    end
end

% Remove missing data from derivatives (represented by NaNs)
if nargout >= 3
    g(~isfinite(g)) = 0;
    if nargout >= 4
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
    subplot(2,4,1); imagesc(f', [0 mxf]); axis xy off; title('f'); colorbar
    subplot(2,4,2); imagesc(mu',[0 mxf]); axis xy off; title('mu'); colorbar
    subplot(2,4,3); imagesc(mu',[0 mxf]); axis xy off; title('nmu'); colorbar        
    subplot(2,4,4); imagesc((f - mu)',[0 mxf]); axis xy off; title('f - mu'); colorbar
    subplot(2,4,5); imagesc(dmu{1}'); axis xy off; title('dmux'); colorbar
    subplot(2,4,6); imagesc(dmu{2}'); axis xy off; title('dmuy'); colorbar
    subplot(2,4,7); imagesc(dmu{3}'); axis xy off; title('dmuz'); colorbar        
else
    subplot(2,4,3); imagesc(mu',[0 mxf]); axis xy off; title('nmu'); colorbar
    subplot(2,4,8); imagesc((f - mu)',[0 mxf]); axis xy off; title('f - nmu'); colorbar
end
drawnow
%==========================================================================    

%==========================================================================
function compare_numerical_derivatives(g,H,f,mu,A,x,tau,Mmu,Mf,speak,method)

% Parameters
dm    = size(x);
B     = get_rigid_basis(dm(3) > 1);               
oq    = A.q;   
Nq    = numel(oq);
vsmu  = sqrt(sum(Mmu(1:3,1:3).^2));

% Modulate gradient w tau
g = -tau*g;

% Function value
R  = spm_dexpm(oq,B);  
M  = Mmu\R*Mf;       
y  = affine_transf(M,x);
fx = meansq_objfun(f,mu,y,A,tau,speak,method);

% Step-size
h = 10.^(-8:1:-1);

%--------------------------------------------------------------------------
% Check gradient and Hessian (diagonal) numerically
% Hessian -> https://math.stackexchange.com/questions/1809060/proof-of-the-second-symmetric-derivative
%--------------------------------------------------------------------------

fprintf('---------------------------\n')
fprintf('Check gradient (com_g vs num_g) and Hessian (com_H vs num_H),\n')
fprintf('for a range of parameter increments (h):\n')

for n=1:Nq % Loop over rigid parameters            
    
    for i=1:numel(h) % Loop over step sizes

        fprintf('         h = %g\n',h(i));
        
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
        
        % Gradient
        num_g = (fxph - fxmh)/(2*h(i));
        
        fprintf('com_g(%i)   = %g\n',n,g(n));    
        fprintf('num_g(%i)   = %g\n',n,num_g);
        
%         % Hessian
%         num_H = (fxph + fxmh - 2*fx)/(h(i)^2);
%     
%         fprintf('com_H(%i,%i) = %g\n',n,n,H(n,n));
%         fprintf('num_H(%i,%i) = %g\n',n,n,num_H);
    end
    
    fprintf('---------------------------\n')
end
%==========================================================================

%==========================================================================
function ok = q_constraint(q,is3d)
mx_tr = 5;
ok    = true;
if is3d
    if q(1) > mx_tr || q(1) < -mx_tr, ok = false; return; end
    if q(2) > mx_tr || q(2) < -mx_tr, ok = false; return; end
    if q(3) > mx_tr || q(3) < -mx_tr, ok = false; return; end
else
    if q(1) > mx_tr || q(1) < -mx_tr, ok = false; return; end
    if q(2) > mx_tr || q(2) < -mx_tr, ok = false; return; end    
end
%==========================================================================

%==========================================================================
function armijo = max_armijo(q,Update,fovmu,is3d)
armijo = 1;
if is3d, imax = 3; else imax = 2; end
for i=1:imax
    armijo = min(armijo, abs((fovmu(i)/2-q(i))/Update(i)));
end
%==========================================================================

function [y,j,d,tol] = cg_im_solver(lhs,rhs,y,nit,tol,verbose)
% Conjugate gradient solver for images
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, y       = zeros(size(rhs),'single'); end
if nargin < 4, nit     = 32; end
if nargin < 5, tol     = 1e-3; end
if nargin < 6, verbose = false; end

% Initilisation  
%--------------------------------------------------------------------------
normRHS = sqrt(sum(rhs(:).*rhs(:))); % Norm of RHS
R       = rhs - lhs(y);              % Residual RHS - LHS(x)
clear RHS

normR = sum(R(:).*R(:));           % R'R
P     = R;                         % Initial conjugate directions P
beta  = 0;                         % Initial search direction for new P

if verbose
    fprintf('0 %g %g\n', tol*normRHS, sqrt(normR));
end

% Run algorithm
%--------------------------------------------------------------------------
j = 1;
while sqrt(normR) > tol*normRHS
    % Calculate conjugate directions P which defines the direction of descent
    %----------------------------------------------------------------------
    P = R + beta*P;

    % Finds the step size of the conj. gradient descent
    %----------------------------------------------------------------------
    AtAP  = lhs(P);
    alpha = normR / sum(P(:).*AtAP(:));
    
    % Perform conj. gradient descent, obtaining updated X and R, using the calculated
    % P and alpha
    %----------------------------------------------------------------------
    y = y + alpha *P; 
    R = R - alpha*AtAP;
    clear alpha AtAP
    
    % Finds the step size for updating P
    %----------------------------------------------------------------------    
    normR = sum(R(:).*R(:));       
    if verbose
        fprintf('%i %g %g\n', j, tol*normRHS, sqrt(normR));
    end
    
    RtRp = normR;    
    beta = normR / RtRp;            
    clear RtRp
    
    % Check if converged
    %----------------------------------------------------------------------
    if j >= nit
        % Finished
        break; 
    end                

    j = j + 1; 
end

d   = sqrt(normR);
tol = tol*normRHS;
%==========================================================================
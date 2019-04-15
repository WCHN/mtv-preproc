function [Nii,ll1,ll3] = update_biasfield(Nii,dat,tau,num_workers,p)
% Update bias field
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Some parameters from options struct
modality      = p.Results.Modality;
speak         = p.Results.Verbose; 
method        = p.Results.Method;
EstimateRigid = p.Results.EstimateRigid;
use_projmat   = ~(strcmpi(method,'denoise') && ~EstimateRigid);
C             = numel(dat); 
if num_workers > 0
    speak = min(speak,1);
end

%--------------------------------------------------------------------------
% Start updating, for each observation
%--------------------------------------------------------------------------

ll1   = zeros(1,C);
ll3   = zeros(1,C);
Nii_b = cell(1,C);
% for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers) % Loop over channels
    
    set_boundary_conditions;
    
    [Nii_b{c},ll1(c),ll3(c)] = update_channel(Nii.x{c},Nii.y(c),Nii.b{c},dat(c),tau{c},use_projmat);    
end % End loop over channels

for c=1:C
    Nii.b{c} = Nii_b{c};
end
clear Nii_b

if speak >= 2
    show_model('bf',dat,Nii,modality);   
end
%==========================================================================

%==========================================================================
function [Nii_b,ll1,ll3] = update_channel(Nii_x,Nii_y,Nii_b,dat,tau,use_projmat) 

N    = numel(Nii_x);
y    = get_nii(Nii_y); 
beta = 1E5;

ll3 = zeros(1,N);
for n=1:N
    
    vx = sqrt(sum(dat.A(n).mat(1:3,1:3).^2));
    
    x        = get_nii(Nii_x(n));
    msk      = get_msk(x);
    bfc      = get_nii(Nii_b(n));
    bf       = exp(bfc);
    bf(~msk) = 1;

    bf        = bf.*y;
    rhs       = tau(n)*bf.*(bf - x);
    lhs       = tau(n)*(bf.^2);
    lhs(~msk) = 0;
    rhs(~msk) = 0;
    clear b
    
    lhs = lhs + spm_field('vel2mom', bfc, [vx 0 0 beta]);
%     lhs = lhs + 1e-7;
    bfc  = spm_field(lhs,rhs,[vx 0 0 beta 2 2]);                      
    clear lhs rhs

    tmp    = spm_field('vel2mom', bfc, [vx 0 0 beta]);
    tmp    = -0.5*bfc(:)'*tmp(:);
    ll3(n) = tmp;

    Nii_b(n) = put_nii(Nii_b(n),bfc);
    clear bc
end

ll1 = get_ll1(use_projmat,true,y,Nii_x,Nii_b,tau,dat);
ll3 = sum(ll3);
%==========================================================================
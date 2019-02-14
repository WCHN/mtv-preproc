function set_boundary_conditions
% Set model boundary conditions
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

spm_field('boundary',1);
pushpull('boundary',1); 
spm_diffeo('boundary',1); 
%==========================================================================   
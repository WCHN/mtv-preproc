function L = laplace_prior(dm,vx)
% Defines a Laplace prior in Fourier space
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

L = @(Y) spm_field('vel2mom', single(Y), [vx 0 1 0]);

% L = zeros(dm);
% if dm(1)>=2
%     tmp        = 1/(vx(1)^2);
%     L(  1,1,1) = L(  1,1,1) + tmp*2;
%     L(  2,1,1) = L(  2,1,1) - tmp;
%     L(end,1,1) = L(end,1,1) - tmp;
% end
% if dm(2)>=2
%     tmp        = 1/(vx(2)^2);
%     L(1,  1,1) = L(1,  1,1) + tmp*2;
%     L(1,  2,1) = L(1,  2,1) - tmp;
%     L(1,end,1) = L(1,end,1) - tmp;
% end
% if dm(3)>=2
%     tmp        = 1/(vx(3)^2);
%     L(1,1,  1) = L(1,1,  1) + tmp*2;
%     L(1,1,  2) = L(1,1,  2) - tmp;
%     L(1,1,end) = L(1,1,end) - tmp;
% end
% L = single(fftn(L));
% L = @(Y) real(ifftn(fftn(Y).*L));
%==========================================================================
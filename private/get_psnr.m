function val = get_psnr(hat,ref)
% Calculates the peak signal-to-noise ratio (PSNR) for the image in array 
% hat, with the image in array ref as the reference.
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

hat = double(hat(:)); 
ref = double(ref(:));

peak = max(hat);
peak = max(peak,max(ref));

RMSE = sqrt(mean((ref - hat).^2)); 

val = 20*log10(peak/RMSE);
%==========================================================================
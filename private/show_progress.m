function show_progress(ll,nii_x,nii_y,dm,nr,nc)
% Show log-likelihood, and noisy and cleaned image
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

figname          = '(SPM) Denoising progress';
f                = findobj('Type', 'Figure', 'Name', figname);
if isempty(f), f = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', f);  

subplot(121)
if numel(ll) >= 3
    plot(ll(3:end),'r-','LineWidth',2); grid on
else
    plot(ll,'r-','LineWidth',2); grid on
end
title('log-likelihood')

C     = numel(nii_x);
img_x = zeros([nr*dm(1) nc*dm(2)],'single');
img_y = zeros([nr*dm(1) nc*dm(2)],'single');
cnt   = 0;
z     = ceil(dm(3)/2);
for c=1:nc
    for r=1:nr        
        rngr = ((r - 1)*dm(1) + 1):r*dm(1);
        rngc = ((c - 1)*dm(2) + 1):c*dm(2);
        
        img_x(rngr,rngc) = get_nii(nii_x(c),z);
        img_y(rngr,rngc) = get_nii(nii_y(c),z);
        
        cnt = cnt + 1;
        if cnt == C, break; end        
    end
end

subplot(122)
imagesc([img_x; max(img_x(:))*ones([20 size(img_x,2)]); img_y]'); axis off image; colormap(gray);
title('Noisy and cleaned')

drawnow;
%==========================================================================
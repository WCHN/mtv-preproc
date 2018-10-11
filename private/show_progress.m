function show_progress(method,ll,nii_x,nii_y,dm,nr,nc,nii_b)
% Show log-likelihood, and input and solved image
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 8
    nii_b = nifti;
end

figname          = '(SPM) MTV progress';
f                = findobj('Type', 'Figure', 'Name', figname);
if isempty(f), f = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', f);  

biascorrect = ~isempty(nii_b(1).dat);

if biascorrect
    subplot(131)
else
    subplot(121)
end
if numel(ll) >= 3
    plot(ll(3:end),'r-','LineWidth',2); grid on
else
    plot(ll,'r-','LineWidth',2); grid on
end
title('log-likelihood')

C = numel(nii_x);

if biascorrect
    subplot(132)
else
    subplot(122)
end
if strcmpi(method,'superres') 
    img_y = zeros([nr*dm(1) nc*dm(2)],'single');
    cnt   = 0;
    z     = ceil(dm(3)/2);
    for c=1:nc
        for r=1:nr        
            rngr = ((r - 1)*dm(1) + 1):r*dm(1);
            rngc = ((c - 1)*dm(2) + 1):c*dm(2);

            img_y(rngr,rngc) = get_nii(nii_y(c),z);

            cnt = cnt + 1;
            if cnt == C, break; end        
        end
    end
    
    imagesc(img_y'); axis off image; colormap(gray);
    title('Input and solved')
else    
    img_x = zeros([nr*dm(1) nc*dm(2)],'single');
    img_y = zeros([nr*dm(1) nc*dm(2)],'single');
    if biascorrect
        img_b = zeros([nr*dm(1) nc*dm(2)],'single');
    end
    cnt   = 0;
    z     = ceil(dm(3)/2);
    for c=1:nc
        for r=1:nr        
            rngr = ((r - 1)*dm(1) + 1):r*dm(1);
            rngc = ((c - 1)*dm(2) + 1):c*dm(2);

            img_x(rngr,rngc) = get_nii(nii_x(c),z);
            img_y(rngr,rngc) = get_nii(nii_y(c),z);
            if biascorrect
                img_b(rngr,rngc) = exp(get_nii(nii_b(c),z));
            end

            cnt = cnt + 1;
            if cnt == C, break; end        
        end
    end
    
    imagesc([img_x; max(img_x(:))*ones([20 size(img_x,2)]); img_y]'); axis off image; colormap(gray); 
    title('Input and solved')
    if biascorrect
        subplot(133)
        imagesc(img_b'); axis off image; colorbar;
        title('Bias field')
    end
end

drawnow;
%==========================================================================
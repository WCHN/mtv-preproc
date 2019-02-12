function show_progress(method,modality,ll,nii_x,nii_y,dm)
% Show log-likelihood, and input and solved image
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

figname          = '(SPM) MTV progress';
f                = findobj('Type', 'Figure', 'Name', figname);
if isempty(f), f = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', f);  

if dm(3) > 1
    nim  = 3;
    dims = [1 2 3; 2 1 1; 3 3 2];
else
    nim  = 1;
    dims = 3;
end

subplot(1,nim+1,1)
plot(ll,'r-','LineWidth',2); grid on
title('log-likelihood')

C  = numel(nii_x);
nr = floor(sqrt(C));
nc = ceil(C/nr);  
    
for d=1:nim
    subplot(1,nim+1,1+d)
%     if strcmpi(method,'superres')
        img_y = zeros([nr*dm(dims(2,d)) nc*dm(dims(3,d))],'single');
        cnt   = 0;
        z     = ceil(dm(dims(1,d))/2);
        for c=1:nc
            for r=1:nr        
                rngr = ((r - 1)*dm(dims(2,d)) + 1):r*dm(dims(2,d));
                rngc = ((c - 1)*dm(dims(3,d)) + 1):c*dm(dims(3,d));

                img_y(rngr,rngc) = fliplr(reshape(get_nii(nii_y(c),z,dims(1,d)), [dm(dims(2,d)) dm(dims(3,d))]));

                cnt = cnt + 1;
                if cnt == C, break; end        
            end
        end

        imagesc(img_y');
%     else    
%         img_x = zeros([nr*dm(dims(2,d)) nc*dm(dims(3,d))],'single');
%         img_y = zeros([nr*dm(dims(2,d)) nc*dm(dims(3,d))],'single');
%         cnt   = 0;
%         z     = ceil(dm(dims(1,d))/2);
%         for c=1:nc
%             for r=1:nr        
%                 rngr = ((r - 1)*dm(dims(2,d)) + 1):r*dm(dims(2,d));
%                 rngc = ((c - 1)*dm(dims(3,d)) + 1):c*dm(dims(3,d));
% 
%                 img_x(rngr,rngc) = fliplr(reshape(get_nii(nii_x{c},z,dims(1,d)), [dm(dims(2,d)) dm(dims(3,d))]));                 
%                 img_y(rngr,rngc) = fliplr(reshape(get_nii(nii_y(c),z,dims(1,d)), [dm(dims(2,d)) dm(dims(3,d))]));
% 
%                 cnt = cnt + 1;
%                 if cnt == C, break; end        
%             end
%         end
% 
%         if strcmpi(modality,'MRI')
%             imagesc([img_x; max(img_x(:))*ones([20 size(img_x,2)]); img_y]');
%         elseif strcmpi(modality,'CT')
%             imagesc([img_x; max(img_x(:))*ones([20 size(img_x,2)]); img_y]',[0 100]);
%         end               
%     end    
    axis off image; colormap(gray); 
    
    if strcmpi(method,'superres')
        title(['Reconstruction (plane ' num2str(d) ')'])
    else
        title(['Denoised (plane ' num2str(d) ')'])
    end
end
    
drawnow;
%==========================================================================
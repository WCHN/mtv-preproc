function show_mtv_scale(mtv_scale)
% Show MTV scaling variable
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

figname          = '(SPM) MTV scaling';
f                = findobj('Type', 'Figure', 'Name', figname);
if isempty(f), f = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', f);  
clf(f)

imagesc3d(mtv_scale); axis off; colormap(parula); colorbar;

title('MTV scaling')
    
drawnow;
%==========================================================================
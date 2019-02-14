function varargout = show_model(varargin)
% Show various parts of the model
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help show_model
    error('Not enough argument. Type ''help spm_prob'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'solution'
        [varargout{1:nargout}] = show_solution(varargin{:});
    case 'll'
        [varargout{1:nargout}] = show_ll(varargin{:});
    case 'mtv'
        [varargout{1:nargout}] = show_mtv(varargin{:});
    case 'rgb'
        [varargout{1:nargout}] = show_rgb(varargin{:});        
    otherwise
        help show_model
        error('Unknown function %s. Type ''help show_model'' for help.', id)
end
%==========================================================================

%==========================================================================
function show_solution(use_projmat,modality,Nii_x,Nii_y)

figname          = '(SPM) MTV solution';
fig              = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);  

C = numel(Nii_x);

if use_projmat
     for c=1:C
        img = single(Nii_y(c).dat(:,:,:));        
        dm  = size(img);
        dm  = [dm 1];
        ix  = floor(dm./2);
        
        if ix(3) > 1
            % 3D
            tmp = squeeze(img(:,:,ix(3)));
            subplot(3,C,(c-1)*C + 1);        
            show_img(tmp,modality);

            tmp = squeeze(img(:,ix(2),:));
            subplot(3,C,(c-1)*C + 2)
            show_img(tmp,modality);

            tmp = squeeze(img(ix(1),:,:));
            subplot(3,C,(c-1)*C + 3)
            show_img(tmp,modality);
        else
            % 2D
            subplot(1,C,c);        
            show_img(img,modality);
        end
    end   
else
    Cc = 2*C;
    for c=1:C
        img0 = single(Nii_x{c}.dat(:,:,:));        
        dm   = size(img0);
        dm   = [dm 1];
        ix0  = floor(dm./2);
        
        img1 = single(Nii_y(c).dat(:,:,:));      
        dm   = size(img1);
        dm   = [dm 1];
        ix1  = floor(dm./2);
        
        if ix0(3) > 1
            % 3D
            tmp = squeeze(img0(:,:,ix0(3)));
            subplot(3,Cc,(c-1)*Cc + 1);
            show_img(tmp,modality);

            tmp = squeeze(img0(:,ix0(2),:));
            subplot(3,Cc,(c-1)*Cc + 2);
            show_img(tmp,modality);

            tmp = squeeze(img0(ix0(1),:,:));
            subplot(3,Cc,(c-1)*Cc + 3);
            show_img(tmp,modality);

            tmp = squeeze(img1(:,:,ix1(3)));
            subplot(3,Cc,(c-1)*Cc + 4);
            show_img(tmp,modality);

            tmp = squeeze(img1(:,ix1(2),:));
            subplot(3,Cc,(c-1)*Cc + 5);
            show_img(tmp,modality);

            tmp = squeeze(img1(ix1(1),:,:));
            subplot(3,Cc,(c-1)*Cc + 6);
            show_img(tmp,modality);
        else
            % 2D
            subplot(2,C,c)
            show_img(img0,modality);
            
            subplot(2,C,C + c)
            show_img(img1,modality);
        end
    end
end

drawnow;
%==========================================================================

%==========================================================================
function show_ll(ll)

figname          = '(SPM) MTV log-likelihood';
fig              = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);  

plot(ll,'r-','LineWidth',2); grid on
title('log-likelihood')

drawnow;
%==========================================================================

%==========================================================================
function show_mtv(mtv)

figname          = '(SPM) MTV prior';
fig              = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);  

imagesc3d(mtv); axis off image xy; colormap(parula); colorbar;
title('Scaling')

drawnow;
%==========================================================================

%==========================================================================
function show_rgb(Nii_y)

figname          = '(SPM) RGB';
fig              = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);  

C    = numel(Nii_y);
dm   = size(Nii_y(1).dat(:,:,:));
is3d = numel(dm) == 3;

Gmag = imgradient3(single(Nii_y(1).dat(:,:,:)));
for c=2:C
    if is3d
        Gmag = cat(4,Gmag,imgradient3(single(Nii_y(c).dat(:,:,:))));
    else
        Gmag = cat(3,Gmag,imgradient3(single(Nii_y(c).dat(:,:,:))));
    end
end

% Make max value <= 1
Gmag = Gmag/max(Gmag(:));
% Gmag = bsxfun(@rdivide,Gmag,sum(Gmag,3));

% Like Gmag' for RGB data
Gmag = permute(Gmag,[2 1 3]);

imagesc3d(Gmag); axis off image xy;
title('RGB')

drawnow;
%==========================================================================

%==========================================================================
% Helper functions
%==========================================================================

%==========================================================================
function show_img(img,modality)
if strcmpi(modality,'CT')
    imagesc(img',[0 100]);
else
    imagesc(img');
end
axis xy off; colormap(gray); 
%==========================================================================
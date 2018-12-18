function Results = GridSearchIXI

% Some folders
dir_out  = './GridSearchIXI';
dir_tmp  = fullfile(dir_out,'TemporaryFiles');
dir_sr   = fullfile(dir_out,'SuperResolved');

if  ~(exist(dir_tmp,'dir') == 7), mkdir(dir_tmp); end 

files = char({'./data/IXI002-Guys-0828-T1.nii', ...
              './data/IXI002-Guys-0828-T2.nii', ...
              './data/IXI002-Guys-0828-PD.nii'});
Nii0  = nifti(files);
C     = numel(Nii0);

% Parameters
VoxelSize           = 1;
DownSampling        = 1/6;
BSplineDegree       = 4;
Verbose             = 1;
WorkersParfor       = Inf;
ZeroMissingValues   = true;
IterMax             = 40;
RegScaleSuperResMRI = 0.01;

% Values for grid-search
IterMax = [20 25 30 35 40 45 50];
ParameterName = 'IterMax';
Parameter     = IterMax;

%% Copy input to new folder and coreg
Nii_ref = coreg_ims(Nii0,dir_tmp);

%% Simulate thick-sliced data
ds = [DownSampling 1 1; 1 DownSampling 1; 1 1 DownSampling];
bs = [BSplineDegree BSplineDegree BSplineDegree  0 0 0];

Nii_in_sr = nifti;
Nii_in_bs = nifti;
for c=1:C
    img0  = Nii_ref(c).dat(:,:,:);
    mat0 = Nii_ref(c).mat;
    dm0  = size(img0);
    D    = diag([ds(c,:) 1]);
    mat  = mat0/D;
    dm   = floor(D(1:3,1:3)*dm0')';
    
    Nii_dat.mat     = mat;
    Nii_dat.dat.dim = dm;
    dat             = init_dat(Nii_dat,mat0,dm0);
    
    % Apply projection matrix
    img = A(img0,dat);
    
%     % Rescale intensities
%     vx0 = sqrt(sum(mat0(1:3,1:3).^2)); 
%     vx  = sqrt(sum(mat(1:3,1:3).^2)); 
%     scl = prod(vx0./vx);
%     img{1} = scl*img{1};
    
    % Save down-sampled image
    [pth,nam,ext] = fileparts(Nii_ref(c).dat.fname);
    nfname        = fullfile(pth,['ds-sr' nam ext]);
    Nii_in_sr(c)     = create_nii(nfname,img{1},mat,[spm_type('float32') spm_platform('bigend')],'downsampled-sr');
    
    % Apply projection matrix
    [x0,y0,z0] = ndgrid(1:dm(1),1:dm(2),1:dm(3));
    
    T = mat0\mat;    

    x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
    y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
    z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

    coeff                       = spm_bsplinc(img0,bs);
    img                         = spm_bsplins(coeff,x1,y1,z1,bs);    
    img(~isfinite(img) | img<0) = 0;
    
    % Save down-sampled image
    [pth,nam,ext] = fileparts(Nii_ref(c).dat.fname);
    nfname        = fullfile(pth,['ds-bs' nam ext]);
    Nii_in_bs(c)  = create_nii(nfname,img,mat,[spm_type('float32') spm_platform('bigend')],'downsampled-bs');
end
clear dat img img0 coeff

%% Reconstruct using bspline
Nii_bs = nifti;
for c=1:C
    img  = Nii_in_bs(c).dat(:,:,:);
    mat0 = Nii_in_bs(c).mat;
    dm0  = size(img);
    D    = diag([ds(c,:) 1]);
    D    = inv(D);
    mat  = mat0/D;
    dm   = floor(D(1:3,1:3)*dm0')';
    
    [x0,y0,z0] = ndgrid(1:dm(1),1:dm(2),1:dm(3));
    
    T = mat0\mat;    

    x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
    y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
    z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

    coeff                       = spm_bsplinc(img,bs);
    img                         = spm_bsplins(coeff,x1,y1,z1,bs);    
    img(~isfinite(img) | img<0) = 0;
    
    [pth,nam,ext] = fileparts(Nii_ref(c).dat.fname);
    nfname        = fullfile(pth,['bs' nam ext]);
    Nii_bs(c)     = create_nii(nfname,img,mat,[spm_type('float32') spm_platform('bigend')],'bs-upsampled');    
        
    V    = spm_vol;
    V(1) = spm_vol(Nii_ref(c).dat.fname);
    V(2) = spm_vol(nfname);
    V    = reslice_ims(V,0,1);
    
    Nii_bs(c) = nifti(V(end).fname);
end
clear img coeff

%% B-spline results
psnr_bs = zeros(1,C); ssim_bs = zeros(1,C); nrmse_bs = zeros(1,C); corr_bs = zeros(1,C);  
for c=1:C    
    ref = Nii_ref(c).dat(:,:,:);
    bs  = Nii_bs(c).dat(:,:,:);
        
    % Mask
    msk       = isfinite(ref) & ref>=0 & isfinite(bs) & bs>=0;    
    ref(~msk) = 0;
    bs(~msk)  = 0;
    
    ssim_bs(1,c)  = ssim(bs,ref); %3d
    corr_bs(1,c)  = corr3(bs,ref); %3d        
    psnr_bs(1,c)  = get_psnr(bs(:),ref(:)); % 1d
    nrmse_bs(1,c) = goodnessOfFit(bs(:),ref(:),'NRMSE'); % 1d
end
clear ref bs

%% Investigate lambda
I      = numel(Parameter);
Nii_sr = cell(1,I);
for i=1:I
    fprintf('.')
    
    OutputDirectory = fullfile(dir_sr,[ParameterName '-' num2str(Parameter(i))]);
    
    if 1
        Nii_sr{i} = spm_mtv_preproc('InputImages',Nii_in_sr,'Method','superres','VoxelSize',VoxelSize,'Verbose',Verbose,'CleanUp',false, ...
                                    'RegScaleSuperResMRI',RegScaleSuperResMRI,'OutputDirectory',OutputDirectory, ...
                                    'CoRegister',false,'WorkersParfor',WorkersParfor, ...
                                    'ZeroMissingValues',ZeroMissingValues,'IterMax',IterMax(i));
                                
        for c=1:C         
            V    = spm_vol;
            V(1) = spm_vol(Nii_ref(c).dat.fname);
            V(2) = spm_vol(Nii_sr{i}(c).dat.fname);
            V    = reslice_ims(V,0,1);

            Nii_sr{i}(c) = nifti(V(end).fname);
        end        
    else
        Nii_sr{i} = nifti;
        for c=1:C
            [~,nam,ext]  = fileparts(Nii_in(c).dat.fname);
            Nii_sr{i}(c) = nifti(fullfile(OutputDirectory,['res_sr_' nam ext]));
        end
    end                       
end
fprintf('\n')

%% Superres results
psnr_lam = zeros(I,C); ssim_lam = zeros(I,C); nrmse_lam = zeros(I,C); corr_lam = zeros(I,C);  
for i=1:I
    for c=1:C    
        ref = Nii_ref(c).dat(:,:,:);
        sr  = Nii_sr{i}(c).dat(:,:,:);

        % Mask
        msk       = isfinite(ref) & ref>=0 & isfinite(sr) & sr>=0;    
        ref(~msk) = 0;
        sr(~msk)  = 0;

        ssim_lam(i,c)  = ssim(sr,ref); % 3d
        corr_lam(i,c)  = corr3(sr,ref); % 3d       
        psnr_lam(i,c)  = get_psnr(sr(:),ref(:)); % 1d
        nrmse_lam(i,c) = goodnessOfFit(sr(:),ref(:),'NRMSE'); % 1d
    end
end
clear ref sr

%% Plot
figure(666)
colors = hsv(3);

subplot(221)
h = semilogx(Parameter,psnr_lam,'-'); hold on
set(h, {'color'}, num2cell(colors, 2));
[mx,ix] = max(psnr_lam,[],1);
h = plot(Parameter(ix),mx,'kx');
h = semilogx(mean(Parameter),psnr_bs,'x'); hold off
set(h, {'color'}, num2cell(colors, 2));
title('PSNR')
xlabel(ParameterName)

subplot(222)
h = semilogx(Parameter,ssim_lam,'-'); hold on
set(h, {'color'}, num2cell(colors, 2));
[mx,ix] = max(ssim_lam,[],1);
h = plot(Parameter(ix),mx,'kx');
h = semilogx(mean(Parameter),ssim_bs,'x'); hold off
set(h, {'color'}, num2cell(colors, 2));
title('SSIM')
xlabel(ParameterName)

subplot(223)
h = semilogx(Parameter,nrmse_lam,'-'); hold on
set(h, {'color'}, num2cell(colors, 2));
[mx,ix] = max(nrmse_lam,[],1);
h = plot(Parameter(ix),mx,'kx');
h = semilogx(mean(Parameter),nrmse_bs,'x'); hold off
set(h, {'color'}, num2cell(colors, 2));
title('NRMSE')
xlabel(ParameterName)

subplot(224)
h = semilogx(Parameter,corr_lam,'-'); hold on
set(h, {'color'}, num2cell(colors, 2));
[mx,ix] = max(corr_lam,[],1);
h = plot(Parameter(ix),mx,'kx');
h = semilogx(mean(Parameter),corr_bs,'x'); hold off
set(h, {'color'}, num2cell(colors, 2));
title('Corr')
xlabel(ParameterName)

% Save figure
saveas(gcf,['mtv-preproc-sr-ixi-' ParameterName '.png']);

% Save results
Results           = struct;
Results.ssim.lam  = ssim_lam;
Results.ssim.bs   = ssim_bs;
Results.nrmse.lam = nrmse_lam;
Results.nrmse.bs  = nrmse_bs;
Results.corr.lam  = corr_lam;
Results.corr.bs   = corr_bs;
Results.psnr.lam  = psnr_lam;
Results.psnr.bs   = psnr_bs;

save(['res' ParameterName '.mat'],'Results');
%==========================================================================

%==========================================================================
function r = corr3(a,b)

% Subtract means 
a = a - mean(a(:));
b = b - mean(b(:));

% Standard Pearson correlation coefficient formula
r = sum(sum(sum(a.*b)))/sqrt(sum(sum(sum(a.*a)))*sum(sum(sum(b.*b))));
%==========================================================================

%==========================================================================
function val = get_psnr(A,B)
peak = max(A);
peak = max(peak,max(B));
A    = double(A); B = double(B);
d    = sum((A(:)-B(:)).^2)/numel(A);
val  = 10*log10(peak*peak/d);
%==========================================================================

%==========================================================================
function [V,ref_ix] = reslice_ims(V,interp,ref_ix)

if nargin<2, interp = 4; end

N = numel(V);
if N==1
    return;
end

if nargin<3
    % Get image with largest volume (for reslicing using this image as
    % reference)
    vol = zeros(N,3);
    for n=1:N
        vx       = spm_misc('vxsize',V(n).mat);
        vol(n,:) = vx.*V(n).dim;
    end
    vol        = prod(vol,2);
    [~,ref_ix] = max(vol);
end

% Set options
matlabbatch{1}.spm.spatial.coreg.write.ref             = {V(ref_ix).fname};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = interp;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap   = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask   = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'res_';        

% Re-slice
ixs       = 1:N;
source_ix = ixs(ixs~=ref_ix);
for n=source_ix
    matlabbatch{1}.spm.spatial.coreg.write.source = {V(n).fname};     

    output_list = spm_jobman('run',matlabbatch);

    delete(V(n).fname);
    V(n) = spm_vol(output_list{1}.rfiles{1});    
end
%==========================================================================
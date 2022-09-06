%%
clear
close all
clc
%% 1. Import and read the image data.
% DAR_img = double(imread('HZ_IMG/HZ_20201111_64Cu-HZ20_new protocol1.tif'));
% Img_name = '20200827-dt-18f-bone/20200827-dt-18f-bone_9';
% Img_name = 'NB_IMG/64Cu_M1TL_S33_MS03';
Img_name = 'P6B1_S69';
DAR_img = double(imread([Img_name,'.tif']));
% DAR_img = double(imread('/Users/penglu/Desktop/p005_RajiSimone_F18-[Phosphor].tif'));
% DAR_img = [DAR_img,DAR_img1(2:end,:)];
% Some imaging instruments implements square-root transformation.
% Therefore, the raw image data should be converted to their real value.
% Convert_scale = 2.3990718e-5;
% DAR_img = (2^16-1-DAR_img).^2*Convert_scale;
% DAR_img = max(DAR_img(:)) - DAR_img;
%% 2. Estimate the mean and standard deviation of the background.
plot_label = 1;
patch_row_size = 10;
patch_col_size = 10;
[bg_mean, bg_std, bg_mask] = estimate_mean_std_gaussian(DAR_img,patch_row_size,patch_col_size,plot_label);

%% 3. Divide the raw image by a pre-set alpha.
bg_std_norm = 10; % this value is normally set between 10 and 15 for our data.
Alpha = bg_std/bg_std_norm;
% Calulate the mean of the background and image divided by Alpha.
bg_mean_norm = double(bg_mean/Alpha);
DAR_img_norm = double(DAR_img/Alpha);

%% 4. Estimate the CNR of the raw image so as to emperically set the regularization parameter for PSF.
sig_mask = estimate_sig_mask(DAR_img_norm); 

Raw_sig = sig_mask.*DAR_img;
Raw_sig(Raw_sig==0)=[];
Raw_bg = bg_mask.*DAR_img;
Raw_bg(Raw_bg==0)=[];

CNR_Raw = (mean(Raw_sig) - mean(Raw_bg))/std(Raw_bg);
%% 5. Intialize the PSF.
kernel_size = 65;
PSF_init = initialize_PSF(kernel_size);
%% 6. Set the parameters in PG-PEM
Maxiter = 100; % max iteration number
% lambda_psf_PGPEM = 4/sqrt(CNR_Raw); % regularization parameter for PSF
lambda_psf_PGPEM = 5; % regularization parameter for PSF
lambda_X_PGPEM = 1.5e-3; % regularization parameter for image
tol_PGPEM = .5e-3; % threshold to stop the iteration
blind_label = 1;
%% 7. Restore the DAR image by PG-PEM.
[X_out_norm,P_out,fun_all]=PG_PEM(DAR_img_norm,PSF_init,bg_std_norm,bg_mean_norm,lambda_psf_PGPEM,lambda_X_PGPEM,Maxiter,tol_PGPEM,blind_label);
X_out = X_out_norm*Alpha; % Multiply the result by alpha. So the result can be compared with the raw image in the same scale.
%% 8. Show the result.
my_viridis; % import colormap
figure;
subplot(211);
imagesc(DAR_img);axis equal;axis off;colormap(viridis_map(70:end,:));title('Raw image','Fontsize',20)
subplot(212);
imagesc(X_out,[0,0.3*max(X_out(:))]);axis equal;axis off;colormap(viridis_map(70:end,:));title('Restored image by PG-PEM','Fontsize',20)

%%
data = uint32(X_out);
% t = Tiff([Img_name,'_restored.tif'],'w');
t = Tiff(['/Users/penglu/Desktop/p005_RajiSimone_F18-[Phosphor]','_restored.tif'],'w');
% Setup tags
% Lots of info here:
% http://www.mathworks.com/help/matlab/ref/tiffclass.html
tagstruct.ImageLength     = size(data,1);
tagstruct.ImageWidth      = size(data,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)
t.write(data);
t.close();
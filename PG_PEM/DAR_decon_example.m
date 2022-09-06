% 1. segment the background part
% 2. Fitting the distribution
% 3. Estimate the image and PSF.
%%
clear
clc
my_viridis;
DAR_img = double(imread('P1B1_S69.tif'));
DAR_img = (65536 - DAR_img).^2*2.3990718e-5;
DAR_img = max(DAR_img(:)) - DAR_img;
[imrow_DAR,imcol_DAR] = size(DAR_img);
%%
patch_row_size = 15;
patch_col_size = 15;
[bg_mean, bg_std, bg_ROI] = estimate_mean_std_gaussian(DAR_img,patch_row_size,patch_col_size,0);
Mag_ratio = bg_std/10;
mean_g = double(bg_mean/Mag_ratio);
std_g = double(bg_std/Mag_ratio);
DAR_img_norm = double(DAR_img/Mag_ratio);

%% estimate normalized std
[Raw_mask,center] = imsegkmeans(single(DAR_img),2);
[~,idx_bg] = min(center);
bg = double(Raw_mask==idx_bg);
Raw_mask = 1-bg;
Raw_mask = bwareaopen(Raw_mask,100);
Raw_mask = ~Raw_mask;
Raw_mask = bwareaopen(Raw_mask,20);
Raw_mask = ~Raw_mask;

bg = bg_ROI;
sig_mask = Raw_mask;
bg_mask = bg; 

Raw_sig = sig_mask.*DAR_img;
Raw_sig(Raw_sig==0)=[];
Raw_bg = bg_mask.*DAR_img;
Raw_bg(Raw_bg==0)=[];

CNR_Raw = (mean(Raw_sig) - mean(Raw_bg))/std(Raw_bg);
lambda = 4/sqrt(CNR_Raw);
%% kernel function
kernel_size = 65 ;
center = [(kernel_size+1)/2,(kernel_size+1)/2];

xx = 1:kernel_size;
yy = 1:kernel_size;
[XX,YY] = meshgrid(xx,yy);
Mag = 1;
kernel = (Mag*(XX - center(1)).^2 + Mag*(YY - center(2)).^2 + 1).^1;
kernel = 1./kernel;
kernel = kernel/sum(kernel(:));

Maxiter = 100;
lambda_psf_proposed = lambda;
lambda_tv = 1e-3; % 1e-3
tol_proposed = 1e-3; % 0.001
%% Implement PG-PEM
tic
[X_out_norm,P_out,fun_all]=PG_PEM(DAR_img_norm,kernel,std_g,mean_g,Maxiter,lambda_psf_proposed,lambda_tv,tol_proposed);
toc
X_out = X_out_norm*Mag_ratio;
%% plot the results
figure;
subplot(1,3,1);
imagesc((DAR_img),[0,max_range]);axis off;colormap(cur_map);axis equal;
colorbar;title('Raw image','Fontsize',15);
subplot(1,3,2);
imagesc((X_out),[0,max_range]);axis off;colormap(cur_map);axis equal;
colorbar;title('PG-PEM','Fontsize',15);
subplot(1,3,3);
imagesc(conv2(X_out,P_out,'same'),[0,max_range]);axis off;colormap(cur_map);axis equal;
colorbar;title('PG-PEM reblurred','Fontsize',15);

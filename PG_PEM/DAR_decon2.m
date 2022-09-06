% 1. segment the background part
% 2. Fitting the distribution
% 3. Estimate the image and PSF.
%%
clear
clc
my_viridis;
DAR_img = double(imread('P1B1_S69.tif'));
% DAR_img = double(imread('mice kidney/20200404-dt-ra-kidney-MS-Xd1.tif'));
% DAR_img = double(imread('NB_20200613_64Cu_M1Tl_MSBend_1.tif'));

% DAR_img = double(imread('20200925-64cucl2-brain-plate2_5left.tif'));

% DAR_img = double(imread('Diane_img/4T1-survival-Ra-amiloride_2.tif'));
% DAR_img = double(imread('20200824-fdg-PC3-tumors_1.tif'));
% DAR_img = double(imread('20200827-dt-18f-bone_2.tif'));
% DAR_img = double(imread('F18-14wks02/F18-14wks02_3.tif'));
% DAR_img = double(imread('20200925-64cucl2-brain-plate2_8.tif'));
% DAR_img = double(imread('211At-EMCD3.tif'));
% DAR_img = double(imread('NB_832019_BioDAge223Ra_A1andA2_HR06_mice3.tif'));

% DAR_img = double(imread('NB_02152019_Autordiography_Patient1_2_corrected1.tif'));
% DAR_img = double(imread('Diane_img/4T1-survival-Ra-amiloride_1.tif'));

% simu_name = 'dar_simu_img_set/simu_img_DAR79_4000_0_1000_9_alpha100';
% DAR_img = double(imread([simu_name,'.tif']));
% DAR_img = double(imread('PSMA/20200522-10um-po210-500min-crop.tif'));
% DAR_img = double(imread('PSMA/223Rarange1_MS01_S19_72h_restored_4.tif'));
% % 
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

Maxiter = 200;
lambda_psf_proposed = lambda;
lambda_psf_SP = lambda;
lambda_psf_RD = 2;
lambda_psf_RL = 6.5;
lambda_tv = 1e-3; % 1e-3
tol_proposed = 1e-3; % 0.001
tol_RD = .3e-3;
tol_RL = 2e-3;
%%
tic
[X_out_norm,P_out,fun_all]=PG_PEM(DAR_img_norm,kernel,std_g,mean_g,Maxiter,lambda_psf_proposed,lambda_tv,tol_proposed);
toc
[X_out_norm_TV,P_out_TV,fun_all_TV]=deblur_gaussian_poisson_EM_new_TV(DAR_img_norm,kernel,std_g,mean_g,Maxiter,lambda_psf_proposed,lambda_tv ,tol_proposed);
[X_out_norm_NONE,P_out_NONE,fun_all_NONE]=deblur_gaussian_poisson_EM_new_TV(DAR_img_norm,kernel,std_g,mean_g,Maxiter,lambda_psf_proposed,0,tol_proposed);
[X_out_norm_SP,P_out_SP,fun_all_SP]=deblur_gaussian_poisson_EM_mod1(DAR_img_norm,kernel,std_g,mean_g,Maxiter,lambda_psf_SP,lambda_tv,tol_proposed);
[X_out_RD,P_out_RD,fun_all_RD]=deblur_RD_EM(DAR_img,kernel,Maxiter,bg_mean,lambda_psf_RD,lambda_tv,tol_RD,'wavelet');
[X_out_blind,P_out_blind,fun_all_blind] = deblur_TV_EM(DAR_img,kernel,Maxiter,bg_mean,lambda_psf_RL,lambda_tv,tol_RL);
%%

X_out = X_out_norm*Mag_ratio;
X_out_NONE = X_out_norm_NONE*Mag_ratio;
X_out_TV = X_out_norm_TV*Mag_ratio;
X_out_SP = X_out_norm_SP*Mag_ratio;

%% plot the results
% figure;
% % subplot(2,3,1);
% plot(1:size(DAR_img,1),DAR_img(:,90),'b',1:size(DAR_img,1),X_out_blind(:,90),'g',1:size(DAR_img,1),X_out(:,90),'r','Linewidth',2);xlim([0,size(DAR_img,1)]);
% % title('Line profile comparison','Fontsize',15);
% legend('raw image','RL','current method');
% xlabel('row number','Fontsize',15);ylabel('Intensity','Fontsize',15);
% figure;
% subplot(1,5,1);
% imagesc((DAR_img),[0,max_range]);axis off;colormap(cur_map);axis equal;
% % colorbar;title('Raw image','Fontsize',15);
% subplot(1,5,2);
% imagesc((X_out_blind),[0,max_range]);axis off;colormap(cur_map);axis equal;
% % colorbar;title('Richard-Lucy algorithm','Fontsize',15);
% subplot(1,5,4);
% imagesc((X_out),[0,max_range]);axis off;colormap(cur_map);axis equal;
% % colorbar;title('Current method','Fontsize',15);
% subplot(1,5,3);
% imagesc(conv2(X_out_blind,P_out_blind,'same'),[0,max_range]);axis off;colormap(cur_map);axis equal;
% % colorbar;title('Richard-Lucy algorithm reblurred','Fontsize',15);
% subplot(1,5,5);
% imagesc(conv2(X_out,P_out,'same'),[0,max_range]);axis off;colormap(cur_map);axis equal;
% colorbar;
% colorbar;title('Current method reblurred','Fontsize',15);

%%
% Raw_total = 3;
% Col_total = 5;
% Fontsize = 18;
% max_range = max([X_out(:),X_out_RD(:),X_out_blind(:)],[],'all');
% cur_map = viridis_map(50:end,:);
% % cur_map = kindl;
% col_num = round(size(DAR_img,2)/2);
% figure;
% row_num = 1:size(DAR_img,1);
% subplot(Raw_total,Col_total,1);
% plot(row_num,DAR_img(row_num,col_num),'b',...
% row_num,X_out_blind(row_num,col_num),'g',...
% row_num,X_out_RD(row_num,col_num),'m',...
% row_num,X_out_SP(row_num,col_num),'k',...
% row_num,X_out(row_num,col_num),'r','Linewidth',2);
% xlim([row_num(1),row_num(end)]);
% % title('Line profile comparison','Fontsize',15);
% legend('raw image','RL','RD','SP','PG-PEM');
% xlabel('row number');ylabel('Intensity');
% set(gca,'Fontsize',Fontsize)
% % figure;
% subplot(Raw_total,Col_total,6);
% imagesc((DAR_img));axis off;colormap(cur_map);
% colorbar;title('Raw image');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,2);
% imagesc((X_out_blind));axis off;colormap(cur_map);
% colorbar;title('RL');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,3);
% imagesc((X_out_RD));axis off;colormap(cur_map);
% colorbar;title('RD');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,4);
% imagesc((X_out_SP));axis off;colormap(cur_map);
% colorbar;title('SP');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,5);
% imagesc((X_out));axis off;colormap(cur_map);
% colorbar;title('PG-PEM');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,7);
% imagesc(conv2(X_out_blind,P_out_blind,'same'));axis off;colormap(cur_map);
% colorbar;title('RL reblurred');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,8);
% imagesc(conv2(X_out_RD,P_out_RD,'same'));axis off;colormap(cur_map);
% colorbar;title('RD reblurred');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,9);
% imagesc(conv2(X_out_SP,P_out_SP,'same'));axis off;colormap(cur_map);
% colorbar;title('SP reblurred');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,10);
% imagesc(conv2(X_out,P_out,'same'));axis off;colormap(cur_map);
% colorbar;title('PG-PEM reblurred');
% set(gca,'Fontsize',Fontsize)
% 
% subplot(Raw_total,Col_total,11);
% imagesc(log(DAR_img));axis off;colormap(cur_map);
% colorbar;title('Raw log-scale');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,12);
% imagesc(log(X_out_blind),[-14,12]);axis off;colormap(cur_map);
% colorbar;title('RL log-scale');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,13);
% imagesc(log(X_out_RD),[-14,12]);axis off;colormap(cur_map);
% colorbar;title('RD log-scale');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,14);
% imagesc(log(X_out_SP),[-14,12]);axis off;colormap(cur_map);
% colorbar;title('SP log-scale');
% set(gca,'Fontsize',Fontsize)
% subplot(Raw_total,Col_total,15);
% imagesc(log(X_out),[-14,12]);axis off;colormap(cur_map);
% colorbar;title('PG-PEM log-scale');
% set(gca,'Fontsize',Fontsize)

%%
Raw_total = 1;
Col_total = 10;
Fontsize = 12;
max_range = max(X_out(:),[],'all');
cur_map = viridis_map(70:end,:);
% thresh = 50;
% cur_map = [plasma_map(end:-1:end-thresh,:);viridis_map(thresh:end,:)];
% cur_map = hot;
% cur_map = plasma_map;
% cur_map = inferno_map;
% cur_map = blackbody_map;
col_num = round(size(DAR_img,2)/2);
% col_num = 200;
% figure;
% row_num = 1:size(DAR_img,1);
% plot(row_num, (DAR_img(row_num,col_num))/max(DAR_img(:)),'-b',...
% row_num,X_out_blind(row_num,col_num)/max(X_out_blind(:)),':g',...
% row_num,X_out_RD(row_num,col_num)/max(X_out_RD(:)),'-.m',...
% row_num,X_out_SP(row_num,col_num)/max(X_out_SP(:)),'--k',...
% row_num,X_out(row_num,col_num)/max(X_out(:)),'-r','Linewidth',1.5);
% xlim([row_num(1),row_num(end)]);
% legend('Raw image','RL','RD','SP','PG-PEM');
% xlabel('Row number');ylabel('Intensity');
% set(gca,'Fontsize',Fontsize)
% grid on

%%
figure;
ha = tight_subplot(Raw_total,Col_total,0.005,0.005,0.005);
axes(ha(1));
imagesc((DAR_img));axis off;colormap(cur_map);
% colorbar;
axes(ha(2));
imagesc((X_out_blind),[0,max_range]);axis off;colormap(cur_map);
axes(ha(3));
imagesc((X_out_RD),[0,max_range]);axis off;colormap(cur_map);
axes(ha(4));
imagesc((X_out_SP),[0,max_range]);axis off;colormap(cur_map);
axes(ha(5));
imagesc((X_out),[0,max_range]);axis off;colormap(cur_map);
% colorbar;
axes(ha(6));
imagesc(log(DAR_img));axis off;colormap(cur_map);
axes(ha(7));
imagesc(log(X_out_blind),[-14,12]);axis off;colormap(cur_map);
axes(ha(8));
imagesc(log(X_out_RD),[-14,12]);axis off;colormap(cur_map);
axes(ha(9));
imagesc(log(X_out_SP),[-14,12]);axis off;colormap(cur_map);
axes(ha(10));
imagesc(log(X_out),[-14,12]);axis off;colormap(cur_map);
%% bg_var
RAW_ROI = DAR_img.*bg_ROI;
RAW_ROI(bg_ROI==0) = [];
std_RAW = std(RAW_ROI);
COV_Raw = std(RAW_ROI)/mean(RAW_ROI);

PGPEM_ROI = X_out.*bg_ROI;
PGPEM_ROI(bg_ROI==0) = [];
std_PGPEM = std(PGPEM_ROI);
COV_PGPEM = std(PGPEM_ROI)/mean(PGPEM_ROI);

PGPEM_TV_ROI = X_out_TV.*bg_ROI;
PGPEM_TV_ROI(bg_ROI==0) = [];
std_PGPEM_TV = std(PGPEM_TV_ROI);
COV_PGPEM_TV = std(PGPEM_TV_ROI)/mean(PGPEM_TV_ROI);

PGPEM_NONE_ROI = X_out_NONE.*bg_ROI;
PGPEM_NONE_ROI(bg_ROI==0) = [];
std_PGPEM_NONE = std(PGPEM_NONE_ROI);
COV_PGPEM_NONE = std(PGPEM_NONE_ROI)/mean(PGPEM_NONE_ROI);

SP_ROI = X_out_SP.*bg_ROI;
SP_ROI(bg_ROI==0) = [];
std_SP = std(SP_ROI);
COV_SP = std(SP_ROI)/mean(SP_ROI);

RD_ROI = X_out_RD.*bg_ROI;
RD_ROI(bg_ROI==0) = [];
std_RD = std(RD_ROI);
COV_RD = std(RD_ROI)/mean(RD_ROI);

RL_ROI = X_out_blind.*bg_ROI;
RL_ROI(bg_ROI==0) = [];
std_RL = std(RL_ROI);
COV_RL = std(RL_ROI)/mean(RL_ROI);
%% CNR

[Raw_mask,center] = imsegkmeans(single(DAR_img),2);
[~,idx_bg] = min(center);
bg = double(Raw_mask==idx_bg);
Raw_mask = 1-bg;
% se = strel('disk',1);
% Raw_mask = imopen(Raw_mask,se);
Raw_mask = bwareaopen(Raw_mask,100);
Raw_mask = ~Raw_mask;
Raw_mask = bwareaopen(Raw_mask,20);
Raw_mask = ~Raw_mask;
% level = multithresh(DAR_img,1);
% Raw_mask = imquantize(DAR_img,level(1));
% Raw_mask = Raw_mask - min(Raw_mask(:));
bg = bg_ROI;
sig_mask = Raw_mask;
bg_mask = bg; 

Raw_sig = sig_mask.*DAR_img;
Raw_sig(Raw_sig==0)=[];
Raw_bg = bg_mask.*DAR_img;
Raw_bg(Raw_bg==0)=[];
CNR_Raw = (mean(Raw_sig) - mean(Raw_bg))/std(Raw_bg);
CR_Raw = mean(Raw_sig)/mean(Raw_bg);
% 
% Blind_mask = double(imsegkmeans(single(X_out_blind),2))-1;
% X_out_blind_CNR = X_out_blind./max(X_out_blind(:));
Blind_mask = Raw_mask;
Blind_sig = sig_mask.*X_out_blind;
Blind_sig(Blind_sig==0)=[];
Blind_bg = bg_mask.*X_out_blind;
Blind_bg(Blind_bg==0)=[];
CNR_Blind = (mean(Blind_sig) - mean(Blind_bg))/std(Blind_bg);
CR_Blind = mean(Blind_sig)/mean(Blind_bg);


RD_sig = sig_mask.*X_out_RD;
RD_sig(RD_sig==0)=[];
RD_bg = bg_mask.*X_out_RD;
RD_bg(RD_bg==0)=[];
CNR_RD = (mean(RD_sig) - mean(RD_bg))/std(RD_bg);
CR_RD = mean(RD_sig)/mean(RD_bg);

Pro_sig = sig_mask.*X_out;
Pro_sig(Pro_sig==0)=[];
Pro_bg = bg_mask.*X_out;
Pro_bg(Pro_bg==0)=[];
CNR_Pro = (mean(Pro_sig) - mean(Pro_bg))/std(Pro_bg);
CR_Pro = mean(Pro_sig)/mean(Pro_bg);

Pro_TV_sig = sig_mask.*X_out_TV;
Pro_TV_sig(Pro_TV_sig==0)=[];
Pro_TV_bg = bg_mask.*X_out_TV;
Pro_TV_bg(Pro_TV_bg==0)=[];
CNR_TV_Pro = (mean(Pro_TV_sig) - mean(Pro_TV_bg))/std(Pro_TV_bg);
CR_TV_Pro = mean(Pro_TV_sig)/mean(Pro_TV_bg);

Pro_NONE_sig = sig_mask.*X_out_NONE;
Pro_NONE_sig(Pro_NONE_sig==0)=[];
Pro_NONE_bg = bg_mask.*X_out_NONE;
Pro_NONE_bg(Pro_NONE_bg==0)=[];
CNR_NONE_Pro = (mean(Pro_NONE_sig) - mean(Pro_NONE_bg))/std(Pro_NONE_bg);
CR_NONE_Pro = mean(Pro_NONE_sig)/mean(Pro_NONE_bg);

SP_sig = sig_mask.*X_out_SP;
SP_sig(SP_sig==0)=[];
SP_bg = bg_mask.*X_out_SP;
SP_bg(SP_bg==0)=[];
CNR_SP = (mean(SP_sig) - mean(SP_bg))/std(SP_bg);
CR_SP = mean(SP_sig)/mean(SP_bg);

%%
% FWHM_RAW = FRC_FWHM_est(DAR_img,0);
% FWHM_PGPEM = FRC_FWHM_est(X_out,0);
% FWHM_SP = FRC_FWHM_est(X_out_SP,0);
% FWHM_RD = FRC_FWHM_est(X_out_RD,0);
% FWHM_RL = FRC_FWHM_est(X_out_blind,0);
%%
% tmp_mask = zeros(size(DAR_img));
% tmp_mask(20:240,20:220) = 1;
% ROI = double(tmp_mask);
% 
% % ROI = double(true_img>1000);
% % ROI = DAR_img > 0;
% save([simu_name,'_restored.mat'],'X_out');
% save([simu_name,'_restored_SP.mat'],'X_out_SP');
% save([simu_name,'_RD_restored.mat'],'X_out_RD');
% save([simu_name,'_RL_restored.mat'],'X_out_blind');
% c0 = sum(true_img.*ROI,'all');
% c1 = sum(X_out.*ROI,'all');
% c2 = sum(X_out_RD.*ROI,'all');
% c3 = sum(X_out_blind.*ROI,'all');
% c0 = 1;
% c1 = 1;
% c2 = 1;
% c3 = 1;
% RMSE1 = (sum((X_out.*ROI*c0/c1 - true_img.*ROI).^2,'all')/sum(ROI,'all')).^.5;
% RMSE2 = (sum((X_out_RD.*ROI*c0/c2 - true_img.*ROI).^2,'all')/sum(ROI,'all')).^.5;
% RMSE3 = (sum((X_out_blind.*ROI*c0/c3 - true_img.*ROI).^2,'all')/sum(ROI,'all')).^.5;
% RMSE4 = (sum((X_out_SP.*ROI*c0/c1 - true_img.*ROI).^2,'all')/sum(ROI,'all')).^.5;

% colorbar;title('Current method reblurred','Fontsize',15);
%% snr
% [psnr_Pro,snr_Pro] = psnr(X_out,true_img,max(true_img(:)));
% [psnr_RL,snr_RL] = psnr(X_out_blind,true_img,max(true_img(:)));
% [psnr_RD,snr_RD] = psnr(X_out_RD,true_img,max(true_img(:)));
% [psnr_SP,snr_SP] = psnr(X_out_SP,true_img,max(true_img(:)));
%% ssim
% ROI = double(true_img>-10);
% ssim_Pro = ssim(X_out.*ROI,true_img.*ROI,'DynamicRange',max(true_img(:)));
% ssim_RL = ssim(X_out_blind.*ROI,true_img.*ROI,'DynamicRange',max(true_img(:)));
% ssim_RD = ssim(X_out_RD.*ROI,true_img.*ROI,'DynamicRange',max(true_img(:)));
% ssim_SP = ssim(X_out_SP.*ROI,true_img.*ROI,'DynamicRange',max(true_img(:)));
%%
% figure;
% subplot(2,3,1);
% plot(fun_all(:,1),'Linewidth',2);xlabel('Iter Nums');ylabel('||X_n_e_w - X_o_l_d||');set(gca,'Fontsize',15)
% subplot(2,3,2);
% plot(fun_all(:,2),'Linewidth',2);xlabel('Iter Nums');ylabel('||h||');set(gca,'Fontsize',15)
% subplot(2,3,3);
% plot(fun_all(:,3),'Linewidth',2);xlabel('Iter Nums');ylabel('Data fidelity');set(gca,'Fontsize',15)
% subplot(2,3,4);
% plot(fun_all(:,4),'Linewidth',2);xlabel('Iter Nums');ylabel('||X||');set(gca,'Fontsize',15)
% subplot(2,3,5);
% plot(fun_all(:,5),'Linewidth',2);xlabel('Iter Nums');ylabel('sum(X*h)');set(gca,'Fontsize',15)
% subplot(2,3,6);
% plot(fun_all(:,6),'Linewidth',2);xlabel('Iter Nums');ylabel('sum(X)');set(gca,'Fontsize',15)
%%
% figure;
% % % row_num = round(size(DAR_img,2)/2)+20;%round(size(DAR_img,2)/2)+20; %176:412;200
% % col_num = round(size(DAR_img,2)/2)+20;
% % % col_num = 110:220;
% % row_num = 44:224;
% % range_idx = row_num*0.042;
% 
% % row_num = round(size(DAR_img,2)/2)+20;%round(size(DAR_img,2)/2)+20; %176:412;200
% col_num = 210:330;
% % col_num = 110:220;
% row_num = 400;
% range_idx = col_num*0.042;
% 
% plot(range_idx, (DAR_img(row_num,col_num))/max(DAR_img(:)),'-b',...
% range_idx,X_out_blind(row_num,col_num)/max(X_out_blind(:)),':g',...
% range_idx,X_out_RD(row_num,col_num)/max(X_out_RD(:)),'-.m',...
% range_idx,X_out_SP(row_num,col_num)/max(X_out_SP(:)),'--k',...
% range_idx,X_out_NONE(row_num,col_num)/max(X_out_NONE(:)),':b',...
% range_idx,X_out_TV(row_num,col_num)/max(X_out_TV(:)),'-.c',...
% range_idx,X_out(row_num,col_num)/max(X_out(:)),'-r','Linewidth',3);
% xlim([range_idx(1),range_idx(end)]);
% h1 = legend({'Raw','RL','RD','SP','NP','TV','PG-PEM',''},'location','eastoutside','NumColumns',2,'Orientation','horizontal');
% set(h1,'box','off')
% xlabel('Distance (mm)');ylabel('Norm. intensity');
% set(gca,'Fontsize',30)
% set(gca,'linewidth',5);
% set(gcf,'color','white');
% set(gca,'TickDir','out')
% box off
% 
% %%
% % 467,56.508064516129,19.999999999997385,135.491935483871
% 
% Markersize = 10;
% Linewidth = 1.5;
% Fontsize = 25;
% % figure;
% % plot(1:7,data_std(1,:),'-s','Linewidth',Linewidth,'Color',[0 0.4470 0.7410],'MarkerSize',Markersize,'MarkerFaceColor',[0 0.4470 0.7410]);
% % hold on;
% % plot(1:7,data_std(2,:),'-s','Linewidth',Linewidth,'Color',[0.8500 0.3250 0.0980],'MarkerSize',Markersize,'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% % hold on;
% % plot(1:7,data_std(3,:),'-s','Linewidth',Linewidth,'Color',[0.9290 0.6940 0.1250],'MarkerSize',Markersize,'MarkerFaceColor',[0.9290 0.6940 0.1250]);
% % hold on;
% % plot(1:7,data_std(4,:),'-s','Linewidth',Linewidth,'Color',[0.4940 0.1840 0.5560],'MarkerSize',Markersize,'MarkerFaceColor',[0.4940 0.1840 0.5560]);
% % hold on;
% % plot(1:7,data_std(5,:),'-s','Linewidth',Linewidth,'Color',[0.4660 0.6740 0.1880],'MarkerSize',Markersize,'MarkerFaceColor',[0.4660 0.6740 0.1880]);
% % hold on;
% % plot(1:7,data_std(6,:),'-s','Linewidth',Linewidth,'Color',[0.3010 0.7450 0.9330],'MarkerSize',Markersize,'MarkerFaceColor',[0.3010 0.7450 0.9330]);
% % hold on;
% % plot(1:7,data_std(7,:),'-s','Linewidth',Linewidth,'Color',[0.6350 0.0780 0.1840],'MarkerSize',Markersize,'MarkerFaceColor',[0.6350 0.0780 0.1840]);
% % hold on;
% % plot(1:7,data_std(8,:),'-s','Linewidth',Linewidth,'Color',[1 0 0],'MarkerSize',Markersize,'MarkerFaceColor',[1 0 0]);
% % hold on;
% % plot(1:7,data_std(9,:),'-s','Linewidth',Linewidth,'Color',[0 1 0],'MarkerSize',Markersize,'MarkerFaceColor',[0 1 0]);
% % hold on;
% % plot(1:7,data_std(10,:),'-s','Linewidth',Linewidth,'Color',[0 0 1],'MarkerSize',Markersize,'MarkerFaceColor',[0 0 1]);
% % 
% % xlim([0,8])
% % xticks([0 1 2 3 4 5 6 7]);
% % xticklabels({'','Raw','RL','RD','SP','NP','TV','PG-PEM'})
% % set(gca,'XTickLabelRotation',30)
% % % xlabel()
% % ylabel('Std')
% % set(gca,'Fontsize',Fontsize)
% % set(gca,'linewidth',5);
% % set(gcf,'color','white');
% % set(gca,'TickDir','out')
% % box off
% 
% data_CNR = [13.27008047	19.94968831	42.59680549	35.12816868	65.72396632	69.55147031	70.18546278
% 12.52721185	18.52382196	46.35926629	35.5900929	65.51728911	69.23348507	69.818717
% 11.85990037	17.68243009	45.54853335	36.22155613	69.09397092	73.6228344	74.74116737
% 9.095168126	10.90530145	27.9903788	24.09404046	47.35893827	50.76730826	51.57293656
% 12.58636839	19.51002364	40.83827517	42.86542932	70.79245511	74.53805263	75.53693141
% 10.66395881	14.40207575	38.06170479	32.54394593	74.49950885	79.9464331	81.01149246
% 8.387595876	9.575464105	24.65626045	28.49017359	49.02324281	51.63987208	51.86187526
% 13.20221242	17.6376215	45.44823157	38.98930327	61.93977651	64.70282775	65.4701141
% 14.18170004	19.93524723	54.82413489	48.72093542	81.58427345	86.06980516	86.87370672
% 11.44060222	13.61013464	41.50607375	32.11774677	53.11442896	56.76712092	57.45046278];
% 
% 
% figure;
% % subplot(321)
% plot(1:7,data_CNR(1,:),'-s','Linewidth',Linewidth,'Color',[0 0.4470 0.7410],'MarkerSize',Markersize,'MarkerFaceColor',[0 0.4470 0.7410]);
% hold on;
% plot(1:7,data_CNR(2,:),'-s','Linewidth',Linewidth,'Color',[0.8500 0.3250 0.0980],'MarkerSize',Markersize,'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% hold on;
% plot(1:7,data_CNR(3,:),'-s','Linewidth',Linewidth,'Color',[0.9290 0.6940 0.1250],'MarkerSize',Markersize,'MarkerFaceColor',[0.9290 0.6940 0.1250]);
% hold on;
% plot(1:7,data_CNR(4,:),'-s','Linewidth',Linewidth,'Color',[0.4940 0.1840 0.5560],'MarkerSize',Markersize,'MarkerFaceColor',[0.4940 0.1840 0.5560]);
% hold on;
% plot(1:7,data_CNR(5,:),'-s','Linewidth',Linewidth,'Color',[0.4660 0.6740 0.1880],'MarkerSize',Markersize,'MarkerFaceColor',[0.4660 0.6740 0.1880]);
% hold on;
% plot(1:7,data_CNR(6,:),'-s','Linewidth',Linewidth,'Color',[0.3010 0.7450 0.9330],'MarkerSize',Markersize,'MarkerFaceColor',[0.3010 0.7450 0.9330]);
% hold on;
% plot(1:7,data_CNR(7,:),'-s','Linewidth',Linewidth,'Color',[0.6350 0.0780 0.1840],'MarkerSize',Markersize,'MarkerFaceColor',[0.6350 0.0780 0.1840]);
% hold on;
% plot(1:7,data_CNR(8,:),'-s','Linewidth',Linewidth,'Color',[1 0 0],'MarkerSize',Markersize,'MarkerFaceColor',[1 0 0]);
% hold on;
% plot(1:7,data_CNR(9,:),'-s','Linewidth',Linewidth,'Color',[0 1 0],'MarkerSize',Markersize,'MarkerFaceColor',[0 1 0]);
% hold on;
% plot(1:7,data_CNR(10,:),'-s','Linewidth',Linewidth,'Color',[0 0 1],'MarkerSize',Markersize,'MarkerFaceColor',[0 0 1]);
% 
% xlim([0,8])
% xticks([0 1 2 3 4 5 6 7]);
% xticklabels({'','Raw','RL','RD','SP','NP','TV','PG-PEM'})
% set(gca,'XTickLabelRotation',45)
% % xlabel()
% ylabel('CNR')
% set(gca,'Fontsize',Fontsize)
% set(gca,'linewidth',5);
% set(gcf,'color','white');
% set(gca,'TickDir','out')
% box off 
% % ax1 = axes('Position',get(gca,'Position'),...
% %            'XAxisLocation','top',...
% %            'YAxisLocation','right',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% % set(ax1,'YTick', []);
% % set(ax1,'XTick', []);
% % box on 
% % set(gca,'linewidth',5);
% 
% data_Res = [0.4863	0.8112	0.7783	0.808	0.8158	0.7952	0.7952
% 0.452	0.8172	0.7943	0.8172	0.8182	0.8184	0.8182
% 0.4571	0.8572	0.8163	0.8555	0.8572	0.8551	0.8587
% 0.4262	0.8087	0.7943	0.803	0.8225	0.8163	0.8081
% 0.483	0.8201	0.8192	0.8188	0.8167	0.8322	0.8309
% 0.452	0.8213	0.7943	0.8542	0.8601	0.8529	0.8592
% 0.4394	0.8095	0.7983	0.8095	0.8251	0.8042	0.8136
% 0.4602	0.8204	0.7879	0.8182	0.8384	0.8182	0.8225
% 0.4123	0.7705	0.6814	0.7755	0.8159	0.7755	0.7755
% 0.4202	0.8151	0.8159	0.8259	0.8531	0.8467	0.8522];
% 
% data_Res = 2*0.042./data_Res;
% 
% % figure;
% % plot(1:7,data_Res(1,:),'-s','Linewidth',Linewidth,'Color',[0 0.4470 0.7410],'MarkerSize',Markersize,'MarkerFaceColor',[0 0.4470 0.7410]);
% % hold on;
% % plot(1:7,data_Res(2,:),'-s','Linewidth',Linewidth,'Color',[0.8500 0.3250 0.0980],'MarkerSize',Markersize,'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% % hold on;
% % plot(1:7,data_Res(3,:),'-s','Linewidth',Linewidth,'Color',[0.9290 0.6940 0.1250],'MarkerSize',Markersize,'MarkerFaceColor',[0.9290 0.6940 0.1250]);
% % hold on;
% % plot(1:7,data_Res(4,:),'-s','Linewidth',Linewidth,'Color',[0.4940 0.1840 0.5560],'MarkerSize',Markersize,'MarkerFaceColor',[0.4940 0.1840 0.5560]);
% % hold on;
% % plot(1:7,data_Res(5,:),'-s','Linewidth',Linewidth,'Color',[0.4660 0.6740 0.1880],'MarkerSize',Markersize,'MarkerFaceColor',[0.4660 0.6740 0.1880]);
% % hold on;
% % plot(1:7,data_Res(6,:),'-s','Linewidth',Linewidth,'Color',[0.3010 0.7450 0.9330],'MarkerSize',Markersize,'MarkerFaceColor',[0.3010 0.7450 0.9330]);
% % hold on;
% % plot(1:7,data_Res(7,:),'-s','Linewidth',Linewidth,'Color',[0.6350 0.0780 0.1840],'MarkerSize',Markersize,'MarkerFaceColor',[0.6350 0.0780 0.1840]);
% % hold on;
% % plot(1:7,data_Res(8,:),'-s','Linewidth',Linewidth,'Color',[1 0 0],'MarkerSize',Markersize,'MarkerFaceColor',[1 0 0]);
% % hold on;
% % plot(1:7,data_Res(9,:),'-s','Linewidth',Linewidth,'Color',[0 1 0],'MarkerSize',Markersize,'MarkerFaceColor',[0 1 0]);
% % hold on;
% % plot(1:7,data_Res(10,:),'-s','Linewidth',Linewidth,'Color',[0 0 1],'MarkerSize',Markersize,'MarkerFaceColor',[0 0 1]);
% % 
% % xlim([0,8])
% % xticks([0 1 2 3 4 5 6 7]);
% % xticklabels({'','Raw','RL','RD','SP','NP','TV','PG-PEM'})
% % set(gca,'XTickLabelRotation',30)
% % % xlabel()
% % ylabel('Resolution')
% % set(gca,'Fontsize',Fontsize)
% % set(gca,'linewidth',5);
% % set(gcf,'color','white');
% % set(gca,'TickDir','out')
% % box off
% 
% data_std = [1318.741616	1008.367085	455.4775288	585.7058797	318.2298516	300.3450242	297.6741441
% 1353.178467	1042.70764	400.5216276	557.2244778	307.7429004	290.8165437	288.4044896
% 1309.64918	1004.394323	371.4358208	503.5057733	268.3793996	251.4782279	247.7239259
% 1329.518944	1253.971941	458.3112664	592.5313507	306.7265883	285.6679557	281.2130665
% 1297.775607	937.9453475	425.7726741	440.6409315	269.6459262	255.9636472	252.6674844
% 1302.266739	1073.750441	391.0693662	497.1074514	219.3589091	204.2333786	201.6279076
% 1308.236455	1261.264209	454.2645357	451.84114	266.1882911	252.1654263	251.0815183
% 1314.433375	1121.897088	423.5390944	543.4537565	346.5252354	331.6792386	328.0585195
% 1335.706708	1089.57712	385.3218526	471.6226712	285.8903185	270.9664337	268.5837462
% 1347.859647	1261.96372	403.48521	573.2100032	351.3802277	328.624973	324.8155698];
% 
% figure;
% yyaxis left
% h2 = boxplot(data_std(:,1),'Widths',0.7,'BoxStyle','outline','Whisker',50,'FullFactors','on');
% % ylim([1296 1360])
% ylabel('Std')
% h22 = findobj(gca,'Tag','Box');
% colors = [0.3010 0.7450 0.9330];
% for j=1:length(h22)
%     patch(get(h22(j),'XData'),get(h22(j),'YData'),colors(j,:),'FaceAlpha',.5);
% end
% set(gca,'ycolor','k');
% yyaxis right
% h3 = boxplot([nan(10,2),data_std(:,2:end)],'Widths',0.7,'BoxStyle','outline','Whisker',5,'FullFactors','on');
% set(h2,'Linewidth',1.5)
% set(h3,'Linewidth',1.5)
% h33 = findobj(gca,'Tag','Box');
% colors = [[0 0 1]
%     [0 0.4470 0.7410]
%     [0.6350 0.0780 0.1840]
%     [0.8500 0.3250 0.0980]
%     [0.9290 0.6940 0.1250]
%     [0.4940 0.1840 0.5560]
%     [0.4660 0.6740 0.1880]
%     [0 0.4470 0.7410]
%     [0.6350 0.0780 0.1840]];
% for j=1:length(h33)-2
%     patch(get(h33(j),'XData'),get(h33(j),'YData'),colors(j,:),'FaceAlpha',.5);
% end
% set(gca,'ycolor','k');
% line([2 2],[0,20000],'LineStyle','--','Color','k','Linewidth',3);
% xlim([0,9])
% xticks([0 1 2 3 4 5 6 7 8]);
% xticklabels({'','Raw','','RL','RD','SP','NP','TV','PG-PEM'})
% set(gca,'XTickLabelRotation',45)
% % legend({'Raw','RL','RD','SP','PG-PEM (no priori)','PG-PEM (TV)','PG-PEM'})
% set(gca,'Fontsize',Fontsize)
% set(gca,'linewidth',5);
% set(gcf,'color','white');
% set(gca,'TickDir','out')
% box off 
% % ax1 = axes('Position',get(gca,'Position'),...
% %            'XAxisLocation','top',...
% %            'YAxisLocation','right',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% % set(ax1,'YTick', []);
% % set(ax1,'XTick', []);
% % box on 
% % set(gca,'linewidth',5);
% 
% figure;
% % subplot(325)
% yyaxis left
% h2 = boxplot(data_Res(:,1),'Widths',0.7,'BoxStyle','outline','Whisker',5,'FullFactors','on');
% ylabel('Resolution (mm)');
% h22 = findobj(gca,'Tag','Box');
% colors = [0.3010 0.7450 0.9330];
% for j=1:length(h22)
%     patch(get(h22(j),'XData'),get(h22(j),'YData'),colors(j,:),'FaceAlpha',.5);
% end
% set(gca,'ycolor','k');
% yyaxis right
% h3 = boxplot([nan(10,2),data_Res(:,2:end)],'Widths',0.7,'BoxStyle','outline','Whisker',5,'FullFactors','on');
% set(h2,'Linewidth',1.5)
% set(h3,'Linewidth',1.5)
% h33 = findobj(gca,'Tag','Box');
% colors = [[0 0 1]
%     [0 0.4470 0.7410]
%     [0.6350 0.0780 0.1840]
%     [0.8500 0.3250 0.0980]
%     [0.9290 0.6940 0.1250]
%     [0.4940 0.1840 0.5560]
%     [0.4660 0.6740 0.1880]
%     [0 0.4470 0.7410]
%     [0.6350 0.0780 0.1840]];
% for j=1:length(h33)
%     patch(get(h33(j),'XData'),get(h33(j),'YData'),colors(j,:),'FaceAlpha',.5);
% end
% set(gca,'ycolor','k');
% line([2 2],[0,20],'LineStyle','--','Color','k','Linewidth',3);
% xlim([0,9])
% xticks([0 1 2 3 4 5 6 7 8]);
% xticklabels({'','Raw','','RL','RD','SP','NP','TV','PG-PEM'})
% set(gca,'XTickLabelRotation',45)
% set(gca,'Fontsize',Fontsize)
% set(gca,'linewidth',5);
% set(gcf,'color','white');
% set(gca,'TickDir','out')
% box off 
% % ax1 = axes('Position',get(gca,'Position'),...
% %            'XAxisLocation','top',...
% %            'YAxisLocation','right',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% % set(ax1,'YTick', []);
% % set(ax1,'XTick', []);
% % box on 
% % set(gca,'linewidth',5);
% % box off
% %%
% data_std2 = [6070.387487	1907.022068	555.4279053	954.2888001	360.6748042	339.4456069	337.9895809
% 6044.603312	2353.923846	599.9028454	902.3709962	240.9672348	224.956208	226.9728471
% 6043.810075	1870.038865	613.7623525	952.9440851	303.4280006	303.099358	302.9524949
% 6094.41863	2189.626747	668.0230895	1049.752067	479.8211685	455.8430782	459.4749636
% 6078.113843	1992.659438	742.6515698	1024.984323	367.5362486	363.4161395	353.9950473
% 6084.386518	1876.600705	447.3187935	953.6146897	299.9915066	299.6766873	293.4835503
% 6052.911165	1898.708362	550.3508801	848.2166085	196.7614478	175.9058783	177.6183402
% 6.101798616854288e+03	1664.836267	615.9398999	1051.498044	559.7911258	540.7506699	545.0195552
% 6.039419470104153e+03	2222.720918	547.6725194	899.7777322	259.2407867	248.7718765	247.9249646
% 6.103375935301427e+03	2150.76449	885.015618	1018.81292	567.1665617	575.1071837	572.7665049];
% 
% 
% % figure;
% % plot(1:7,data_std2(1,:),'-s','Linewidth',Linewidth,'Color',[0 0.4470 0.7410],'MarkerSize',Markersize,'MarkerFaceColor',[0 0.4470 0.7410]);
% % hold on;
% % plot(1:7,data_std2(2,:),'-s','Linewidth',Linewidth,'Color',[0.8500 0.3250 0.0980],'MarkerSize',Markersize,'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% % hold on;
% % plot(1:7,data_std2(3,:),'-s','Linewidth',Linewidth,'Color',[0.9290 0.6940 0.1250],'MarkerSize',Markersize,'MarkerFaceColor',[0.9290 0.6940 0.1250]);
% % hold on;
% % plot(1:7,data_std2(4,:),'-s','Linewidth',Linewidth,'Color',[0.4940 0.1840 0.5560],'MarkerSize',Markersize,'MarkerFaceColor',[0.4940 0.1840 0.5560]);
% % hold on;
% % plot(1:7,data_std2(5,:),'-s','Linewidth',Linewidth,'Color',[0.4660 0.6740 0.1880],'MarkerSize',Markersize,'MarkerFaceColor',[0.4660 0.6740 0.1880]);
% % hold on;
% % plot(1:7,data_std2(6,:),'-s','Linewidth',Linewidth,'Color',[0.3010 0.7450 0.9330],'MarkerSize',Markersize,'MarkerFaceColor',[0.3010 0.7450 0.9330]);
% % hold on;
% % plot(1:7,data_std2(7,:),'-s','Linewidth',Linewidth,'Color',[0.6350 0.0780 0.1840],'MarkerSize',Markersize,'MarkerFaceColor',[0.6350 0.0780 0.1840]);
% % hold on;
% % plot(1:7,data_std2(8,:),'-s','Linewidth',Linewidth,'Color',[1 0 0],'MarkerSize',Markersize,'MarkerFaceColor',[1 0 0]);
% % hold on;
% % plot(1:7,data_std2(9,:),'-s','Linewidth',Linewidth,'Color',[0 1 0],'MarkerSize',Markersize,'MarkerFaceColor',[0 1 0]);
% % hold on;
% % plot(1:7,data_std2(10,:),'-s','Linewidth',Linewidth,'Color',[0 0 1],'MarkerSize',Markersize,'MarkerFaceColor',[0 0 1]);
% % 
% % xlim([0,8])
% % xticks([0 1 2 3 4 5 6 7]);
% % xticklabels({'','Raw','RL','RD','SP','NP','TV','PG-PEM'})
% % set(gca,'XTickLabelRotation',30)
% % % xlabel()
% % ylabel('Std')
% % set(gca,'Fontsize',Fontsize)
% % set(gca,'linewidth',5);
% % set(gcf,'color','white');
% % set(gca,'TickDir','out')
% % box off
% % % figure;
% % % boxplot(data_std2,'Widths',0.5);
% 
% data_CNR2 = [6.86617335	24.11587009	75.52771705	47.24309242	127.4416156	135.3248001	136.0234854
% 4.200239409	12.46123128	41.97156071	31.58245774	119.712829	128.012447	127.0314868
% 5.635322149	22.97183039	59.5440861	42.14748508	135.8792681	135.2859539	135.5698756
% 6.696107143	22.50713095	62.95813595	45.14520867	100.9864977	106.1201956	105.3630482
% 7.853961776	28.04790479	66.93790701	51.23143238	146.0168043	147.1495018	151.3901781
% 7.619990946	27.33730147	107.6384692	53.46081186	173.3390186	172.9093742	176.8502762
% 4.429268848	15.09213492	48.39729036	32.6687675	145.9459529	163.5487233	162.0368442
% 6.840276378	30.69516714	72.25235398	44.89445567	86.34525459	89.18452669	88.54934768
% 4.663618462	14.70264839	51.73568431	34.97512642	124.3940815	129.178378	129.8560751
% 5.603833195	20.50491246	43.13421193	39.76363588	72.96295432	71.44456095	71.86615131
%     ];
% 
% 
% figure;
% plot(1:7,data_CNR2(1,:),'-s','Linewidth',Linewidth,'Color',[0 0.4470 0.7410],'MarkerSize',Markersize,'MarkerFaceColor',[0 0.4470 0.7410]);
% hold on;
% plot(1:7,data_CNR2(2,:),'-s','Linewidth',Linewidth,'Color',[0.8500 0.3250 0.0980],'MarkerSize',Markersize,'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% hold on;
% plot(1:7,data_CNR2(3,:),'-s','Linewidth',Linewidth,'Color',[0.9290 0.6940 0.1250],'MarkerSize',Markersize,'MarkerFaceColor',[0.9290 0.6940 0.1250]);
% hold on;
% plot(1:7,data_CNR2(4,:),'-s','Linewidth',Linewidth,'Color',[0.4940 0.1840 0.5560],'MarkerSize',Markersize,'MarkerFaceColor',[0.4940 0.1840 0.5560]);
% hold on;
% plot(1:7,data_CNR2(5,:),'-s','Linewidth',Linewidth,'Color',[0.4660 0.6740 0.1880],'MarkerSize',Markersize,'MarkerFaceColor',[0.4660 0.6740 0.1880]);
% hold on;
% plot(1:7,data_CNR2(6,:),'-s','Linewidth',Linewidth,'Color',[0.3010 0.7450 0.9330],'MarkerSize',Markersize,'MarkerFaceColor',[0.3010 0.7450 0.9330]);
% hold on;
% plot(1:7,data_CNR2(7,:),'-s','Linewidth',Linewidth,'Color',[0.6350 0.0780 0.1840],'MarkerSize',Markersize,'MarkerFaceColor',[0.6350 0.0780 0.1840]);
% hold on;
% plot(1:7,data_CNR2(8,:),'-s','Linewidth',Linewidth,'Color',[1 0 0],'MarkerSize',Markersize,'MarkerFaceColor',[1 0 0]);
% hold on;
% plot(1:7,data_CNR2(9,:),'-s','Linewidth',Linewidth,'Color',[0 1 0],'MarkerSize',Markersize,'MarkerFaceColor',[0 1 0]);
% hold on;
% plot(1:7,data_CNR2(10,:),'-s','Linewidth',Linewidth,'Color',[0 0 1],'MarkerSize',Markersize,'MarkerFaceColor',[0 0 1]);
% 
% xlim([0,8])
% xticks([0 1 2 3 4 5 6 7]);
% xticklabels({'','Raw','RL','RD','SP','NP','TV','PG-PEM'})
% set(gca,'XTickLabelRotation',45)
% % xlabel()
% % ylabel('CNR')
% set(gca,'Fontsize',Fontsize)
% set(gca,'linewidth',5);
% set(gcf,'color','white');
% set(gca,'TickDir','out')
% box off 
% % ax1 = axes('Position',get(gca,'Position'),...
% %            'XAxisLocation','top',...
% %            'YAxisLocation','right',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% % set(ax1,'YTick', []);
% % set(ax1,'XTick', []);
% % box on 
% % set(gca,'linewidth',5);
% 
% figure;
% h2 = boxplot(data_std2,'Widths',0.7,'BoxStyle','outline','Whisker',5,'FullFactors','on');
% breakyaxis([2500,5000])
% % ylabel('Std');
% h22 = findobj(gca,'Tag','Box');
% set(h2,'Linewidth',1.5)
% colors = [
%     [0 0.4470 0.7410]
%     [0.6350 0.0780 0.1840]
%     [0.8500 0.3250 0.0980]
%     [0.9290 0.6940 0.1250]
%     [0.4940 0.1840 0.5560]
%     [0.4660 0.6740 0.1880]
%     [0 0.4470 0.7410]
%     [0.6350 0.0780 0.1840]];
% for j=1:length(h22)
%     patch(get(h22(j),'XData'),get(h22(j),'YData'),colors(j,:),'FaceAlpha',.5);
% end
% xlim([0,8])
% xticks([0 1 2 3 4 5 6 7]);
% xticklabels({'','Raw','RL','RD','SP','NP','TV','PG-PEM'})
% set(gca,'XTickLabelRotation',45)
% % legend({'Raw','RL','RD','SP','PG-PEM (no priori)','PG-PEM (TV)','PG-PEM'})
% set(gca,'Fontsize',Fontsize)
% set(gca,'linewidth',5);
% set(gcf,'color','white');
% set(gca,'TickDir','out')
% box off 
% % ax1 = axes('Position',get(gca,'Position'),...
% %            'XAxisLocation','top',...
% %            'YAxisLocation','right',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% % set(ax1,'YTick', []);
% % set(ax1,'XTick', []);
% % box on 
% % set(gca,'linewidth',5);
% 
% 
% data_Res2 = [0.03343	0.348	0.3429	0.3665	0.3479	0.3561	0.3768
% 0.03343	0.3565	0.3552	0.3566	0.3528	0.3552	0.3559
% 0.05335	0.3521	0.3426	0.3507	0.3376	0.3384	0.3532
% 0.06904	0.3521	0.3436	0.3573	0.3236	0.3573	0.3569
% 0.03394	0.354	0.3158	0.3531	0.3536	0.3557	0.3564
% 0.04941	0.3453	0.3407	0.3454	0.3455	0.3455	0.3489
% 0.02511	0.3597	0.353	0.3528	0.3316	0.3348	0.3606
% 0.06552	0.3632	0.3146	0.348	0.3184	0.3192	0.3599
% 0.03343	0.354	0.3446	0.3549	0.3516	0.3479	0.3557
% 0.06929	0.366	0.3465	0.3429	0.35	0.3408	0.3506];
% 
% data_Res2 = 2*0.042./data_Res2;
% 
% figure;
% yyaxis left
% h2 = boxplot(data_Res2(:,1),'Widths',0.7,'BoxStyle','outline','Whisker',5,'FullFactors','on');
% % ylabel('Resolution (mm)');
% h22 = findobj(gca,'Tag','Box');
% colors = [0.3010 0.7450 0.9330];
% for j=1:length(h22)
%     patch(get(h22(j),'XData'),get(h22(j),'YData'),colors(j,:),'FaceAlpha',.5);
% end
% set(gca,'ycolor','k');
% yyaxis right
% h3 = boxplot([nan(10,2),data_Res2(:,2:end)],'Widths',0.7,'BoxStyle','outline','Whisker',5,'FullFactors','on');
% set(h2,'Linewidth',1.5)
% set(h3,'Linewidth',1.5)
% h33 = findobj(gca,'Tag','Box');
% colors = [[0 0 1]
%     [0 0.4470 0.7410]
%     [0.6350 0.0780 0.1840]
%     [0.8500 0.3250 0.0980]
%     [0.9290 0.6940 0.1250]
%     [0.4940 0.1840 0.5560]
%     [0.4660 0.6740 0.1880]
%     [0 0.4470 0.7410]
%     [0.6350 0.0780 0.1840]];
% for j=1:length(h33)
%     patch(get(h33(j),'XData'),get(h33(j),'YData'),colors(j,:),'FaceAlpha',.5);
% end
% set(gca,'ycolor','k');
% line([2 2],[0,20],'LineStyle','--','Color','k','Linewidth',3);
% xlim([0,9])
% xticks([0 1 2 3 4 5 6 7 8]);
% xticklabels({'','Raw','','RL','RD','SP','NP','TV','PG-PEM'})
% set(gca,'XTickLabelRotation',45)
% set(gca,'Fontsize',Fontsize)
% set(gca,'linewidth',5);
% set(gcf,'color','white');
% set(gca,'TickDir','out')
% box off 
% % ax1 = axes('Position',get(gca,'Position'),...
% %            'XAxisLocation','top',...
% %            'YAxisLocation','right',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% % set(ax1,'YTick', []);
% % set(ax1,'XTick', []);
% % box on 
% % set(gca,'linewidth',5);

%%
% figure;plot(1:length(fun_all(:,3)),fun_all(:,3),'Linewidth',2,'Color',[0.9290 0.6940 0.1250]);
% ylabel('Log-likelihood')
% set(gca,'Fontsize',25)
% set(gca,'linewidth',5);
% set(gcf,'color','white');
% set(gca,'TickDir','out')
% box off
% 
% figure;plot(2:length(fun_all(:,7)),fun_all(2:end,7),'Linewidth',2,'Color',[0.4940 0.1840 0.5560]);
% ylabel('||X_{new}-X_{old}||/||X_{old}||');
% set(gca,'Fontsize',25)
% set(gca,'linewidth',5);
% set(gcf,'color','white');
% set(gca,'TickDir','out')
% box off
%%
% F_blind = log(abs(fftshift(fft2(X_out_blind))));
% F_RD = log(abs(fftshift(fft2(X_out_RD))));
% F_SP = log(abs(fftshift(fft2(X_out_SP))));
% F_NONE = log(abs(fftshift(fft2(X_out_NONE))));
% F_TV = log(abs(fftshift(fft2(X_out_TV))));
% F_HF = log(abs(fftshift(fft2(X_out))));
% min_range = min([F_blind,F_RD,F_SP,F_NONE,F_TV,F_HF],[],'all')-1;
% max_range = max([F_blind,F_RD,F_SP,F_NONE,F_TV,F_HF],[],'all')+1;
% cur_map = gray;
% 
% figure;
% imagesc(linspace(-1,1,size(F_blind,1)),linspace(-1,1,size(F_blind,2)),log(abs(fftshift(fft2(DAR_img)))));colormap(cur_map);axis equal;axis off
% 
% figure;
% subplot(231)
% imagesc(linspace(-1,1,size(F_blind,1)),linspace(-1,1,size(F_blind,2)),F_blind,[min_range,max_range]);colormap(cur_map);axis equal;axis off
% 
% subplot(232)
% imagesc(linspace(-1,1,size(F_RD,1)),linspace(-1,1,size(F_RD,2)),F_RD,[min_range,max_range]);colormap(cur_map);axis equal;axis off
% 
% subplot(233)
% imagesc(linspace(-1,1,size(F_SP,1)),linspace(-1,1,size(F_SP,2)),F_SP,[min_range,max_range]);colormap(cur_map);axis equal;axis off
% 
% subplot(234)
% imagesc(linspace(-1,1,size(F_NONE,1)),linspace(-1,1,size(F_NONE,2)),F_NONE,[min_range,max_range]);colormap(cur_map);axis equal;axis off
% 
% subplot(235)
% imagesc(linspace(-1,1,size(F_TV,1)),linspace(-1,1,size(F_TV,2)),F_TV,[min_range,max_range]);colormap(cur_map);axis equal;axis off
% 
% subplot(236)
% imagesc(linspace(-1,1,size(F_HF,1)),linspace(-1,1,size(F_HF,2)),F_HF,[min_range,max_range]);colormap(cur_map);axis equal;axis off
function [mean_gaussian, std_gaussian, ROI] = estimate_mean_std_gaussian(Input_img,patch_row_size,patch_col_size,plot_label)
%% patch-based estimation
patch_row = floor(size(Input_img,1)/patch_row_size);
patch_col = floor(size(Input_img,2)/patch_col_size);

mean_list = [];
std_list = [];
skewness_list = [];
kurtosis_list = [];

for ii = 1:patch_row
    for jj = 1:patch_col
        Input_img_patch = Input_img((ii-1)*patch_row_size+1:ii*patch_row_size,(jj-1)*patch_col_size+1:jj*patch_col_size);
        mean_list(end+1) = mean(Input_img_patch(:));
        std_list(end+1) = std(Input_img_patch(:));
        skewness_list(end+1) = skewness(Input_img_patch(:));
        kurtosis_list(end+1) = kurtosis(Input_img_patch(:));
    end
end
%% cluster
mean_list_norm = zscore(mean_list);
std_list_norm = zscore(std_list);
skewness_list_norm = zscore(skewness_list);
kurtosis_list_norm = zscore(kurtosis_list);
%%
[XX_mean,YY_mean] = meshgrid(mean_list_norm,mean_list_norm);
[XX_std,YY_std] = meshgrid(std_list_norm,std_list_norm);
[XX_skewness,YY_skewness] = meshgrid(skewness_list_norm,skewness_list_norm);
[XX_kurtosis,YY_kurtosis] = meshgrid(kurtosis_list_norm,kurtosis_list_norm);
ZZ = sqrt((XX_mean-YY_mean).^2+(XX_std-YY_std).^2+(XX_skewness-YY_skewness).^2+(XX_kurtosis-YY_kurtosis).^2);
ZZ = sort(ZZ,1);
ZZZ = sort(ZZ(8,:));
%%
label = dbscan([mean_list_norm;std_list_norm;skewness_list_norm;kurtosis_list_norm]',ZZZ(floor(size(ZZ,2)*0.25)),8);
label_set = label == 1;

% figure;
% my_viridis;
% max_range = max(Input_img(:));
% ha = tight_subplot(patch_row,patch_col,0.005,0.005,0.005);
% for ii = 1:patch_row
%     for jj = 1:patch_col
% %         if label_set((ii-1)*patch_col+jj) == 1
%             axes(ha((ii-1)*patch_col+jj));
%             Input_img_patch = Input_img((ii-1)*patch_row_size+1:ii*patch_row_size,(jj-1)*patch_col_size+1:jj*patch_col_size);
%             imagesc(Input_img_patch,[0,max_range]);axis off;colormap(viridis_map(70:end,:));
% %         end 
%     end
% end
%%
ROI = zeros(size(Input_img));
for ii = 1:patch_row
    for jj = 1:patch_col
        if label_set((ii-1)*patch_col+jj) == 1
            ROI((ii-1)*patch_row_size+1:ii*patch_row_size,(jj-1)*patch_col_size+1:jj*patch_col_size) = 1;
        end 
    end
end

%%
Input_img_collect = [];
for kk = 1:length(label_set)
    if label_set(kk) == 1
        row_idx = floor((kk-1)/patch_col)+1;
        col_idx = mod(kk-1,patch_col)+1;
        Img_patch = Input_img((row_idx-1)*patch_row_size+1:row_idx*patch_row_size,(col_idx-1)*patch_col_size+1:col_idx*patch_col_size);
        Input_img_collect = cat(3,Input_img_collect,Img_patch);
    end 
end
if plot_label == 1
    mean_list_norm_bg = mean_list_norm(label_set==1);
    std_list_norm_bg = std_list_norm(label_set==1);
    skewness_list_norm_bg = skewness_list_norm(label_set==1);
    mean_list_norm_sig = mean_list_norm(label_set==0);
    std_list_norm_sig = std_list_norm(label_set==0);
    skewness_list_norm_sig = skewness_list_norm(label_set==0);

    figure;
    subplot(1,2,1);
    scatter3(mean_list_norm_bg,std_list_norm_bg,skewness_list_norm_bg,100,'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75]);
    hold on;
    scatter3(mean_list_norm_sig,std_list_norm_sig,skewness_list_norm_sig,100,'MarkerEdgeColor','r',...
            'MarkerFaceColor','b');
    xlabel('Normalized mean');ylabel('Normalized std');zlabel('Normalized skewness');legend({'Background','Signal'})
    set(gca,'Fontsize',30);box on;ax = gca;
    ax.BoxStyle = 'full';colormap(jet)
    title('Statistics distribution of patches');
    set(gca,'Fontsize',15);
    subplot(1,2,2);
    h = histogram(Input_img_collect(:),200);xlabel('Intensity');ylabel('Counts');set(gca,'Fontsize',15)
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'b';
%     title('Histogram of the combined patches');
    set(gca,'Fontsize',15)
end


%% Trucated MLE
dat_normal = Input_img_collect(:);
dat_normal(dat_normal==0) = [];
x_min = min(dat_normal);
[~, phat]  = fitdist_ntrunc(dat_normal, [x_min, Inf]);
mean_gaussian = phat(1);
std_gaussian = phat(2);

% y = 0:1:8000;
% f = exp(-(y-mean_gaussian).^2./(2*std_gaussian^2))./(std_gaussian*sqrt(2*pi));

% figure;
% h = histogram(Input_img(:),121,'Normalization','pdf','LineWidth',1);
% h.FaceColor = [0 1 0];
% h.EdgeColor = 'k';
% 
% box off;
% axis off;
% 
% figure;
% h = histogram(Input_img_collect(:),31,'Normalization','pdf','LineWidth',2);
% % h.FaceColor = [0 1 0];
% h.EdgeColor = 'b';
% hold on;
% plot(y, f, 'Color', [1, 0, 0], 'Linewidth', 3)
% box off;
% axis off;
end
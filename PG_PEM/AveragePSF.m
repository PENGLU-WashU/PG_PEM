function P_out = AveragePSF(P_in)
    pox = (size(P_in,1)+1)/2;
    poy = (size(P_in,2)+1)/2;
    [POX, POY] = meshgrid(1:size(P_in,1),1:size(P_in,2));
    POD = (POX - pox).^2 + (POY - poy).^2;
    POD_list = unique(POD);
    p_tmp = P_in(:);
    for kk = 1:length(POD_list)
        pos_all = find(POD == POD_list(kk));
        p_tmp(pos_all) = (mean((p_tmp(pos_all)).^1,'all')).^1;
    end
    P_out = reshape(p_tmp,size(P_in));
end
% [X_out,P_out,fun_all]=PG_PEM(X_raw,P,sigma,c,lambda_psf,lambda_X,Maxiter,tol,blind_label)
% ---------------------------------------
%
% Blind restore DAR images using Poisson-Gaussian Penalized Maximum
% Likelihood Expectation Maximization Algorithm (PG-PEM).
%
% Inputs:
%  X_raw        Raw image to be restored
%  P           	Initial PSF
%  sigma		Standard deviation of the background
%  c            Mean of the background
%  lambda_psf   Regularization parameter of the PSF, normally 0 to 4
%  lambda_X     Regularization parameter of the image, normally around 1e-3
%  Maxiter      Maximum iteration number
%  tol          Tolerance of the iterations
%  blind_label  if blind_label = 1, blind restore the raw image
%               else restore the raw image with known psf
%
% Outputs:
%  X_out        Restored image
%  P_out		Estimated PSF
%  fun_all		Saved ||X_new-X_old||/||X_old|| in each iterations
%
% ---------------------------------------

function [X_out,P_out,fun_all]=PG_PEM(X_raw,P,sigma,c,lambda_psf,lambda_X,Maxiter,tol,blind_label)
if nargin < 9
    blind_label = 1;
end
if nargin < 8 
    tol = 0.8e-3; 
end
if nargin < 7
    Maxiter = 150;  
end
if nargin < 6 
    lambda_X = 1e-3; 
end
%% parameters initialization
sigma2 = sigma^2;
X_new = X_raw;
P_new = P;
q_max_old = -1;
Xraw_list = [];
[n_rows,n_cols] = size(X_raw);
fun_all = NaN(Maxiter,1);

%% Implementation
for ii = 1:Maxiter 
%% Update PSF and image
    X_old = X_new;
    P_old = P_new;
%% E step    
    switch ii
        case 1 
            U = X_raw;
        otherwise
            U = conv2(X_old,P_old,'same');
    end
     
    qstar = real(n_star(U,X_raw,c,sigma2)); 
    q_ub = round(qstar + 5*sigma); 
    q_max = max(q_ub(:));
    
    if q_max_old < q_max
        Xraw_list_new = NaN(n_rows,n_cols,q_max-q_max_old);
        for qs = (q_max_old+1):q_max
            Xraw_list_new(:,:,qs-q_max_old) = exp(-(X_raw - c - qs).^2/(2*sigma2)); 
        end
        Xraw_list = cat(2,Xraw_list,reshape(Xraw_list_new,[n_rows*n_cols,q_max-q_max_old]));
        q_max_old = q_max;
    end

    E_XH = Calculate_Estep(U,Xraw_list,q_max);
    
%% M step: Update the image
    H = psf2otf(P_new,size(X_raw));
    scale = real(ifftn(conj(H).*fftn(ones(size(H))))) + sqrt(eps);
    tmp = E_XH./(ifftn(fftn(X_old).*H)+eps);
    X_new = X_old.*real(ifftn(conj(H).*fftn(tmp)))./scale;  
    HF_Reg = HessianFrobenious(X_old);
    X_new = X_new./(1+lambda_X*HF_Reg);
    X_new = max(X_new,eps);
    
%% M step: update the PSF
    if blind_label == 1
        H = fftn(X_new);
        scale = otf2psf(conj(H).*fftn(ones(size(H))),size(P_old)) + sqrt(eps);
        tmp = E_XH./(ifftn(H.*psf2otf(P_old,size(X_raw)))+eps);
        P_new = P_old.*real((otf2psf(conj(H).*fftn(tmp),size(P_old)))./scale);
        P_new = P_new./(1+lambda_psf*P_old);
        P_new = AveragePSF(P_new);
        P_new = P_new/sum(P_new(:));
        if ~isfinite(P_new)
            disp('The psf is infinite!');
            break;
        end
    end
%% Check whether the iteration can be stopped
    if ii>1
        fun_all(ii,1) = norm(X_new - X_old)/norm(X_old);
        if fun_all(ii,1) < tol
            break;
        end
    end
    
end

fun_all((ii+1):end,:) = [];
X_out = X_new;
P_out = P_new;
end

function res = n_star(a,b,c,sigma2)
    res = sigma2/1.^2*lambertw(1/sigma2.*a.*exp(1/sigma2*(b-c)));
end

function E_XH = Calculate_Estep(U,Xraw_list,q_max)
    divider = max(exp(-U/2),4.9407e-323);
    omega = 0;
    [n_rows,n_cols] = size(U);
    
    eta = reshape(Xraw_list(:,1),[n_rows,n_cols]).*U.*divider;
    karma = divider;
    for qs = 1:q_max
        switch qs
            case 700
                karma = karma.*U/qs;
                karmast700 = karma.*divider;
                karma = karmast700.*U/qs;
                kesei = karma.*reshape(Xraw_list(:,qs+1),[n_rows,n_cols]);
                omega = omega.*divider+kesei*qs;
                eta = eta.*divider+kesei;
            otherwise
                karma = karma.*U/qs;
                kesei = karma.*reshape(Xraw_list(:,qs+1),[n_rows,n_cols]);
                omega = omega+kesei*qs;
                eta = eta+kesei;
        end
    end
    
    E_XH = omega./(eta+eps);
end
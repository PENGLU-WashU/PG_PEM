function [X_out,P_out,fun_all]=PG_PEM(Bobs,P,sigma,c,Maxiter,lambda_psf,lambda_HF,tol)
%% parameters initialization
sigma2 = sigma^2;
X_new = Bobs;
fun_all = NaN(Maxiter,3);
P_new = P;
    
%% Implementation
for ii = 1:Maxiter
    Sst = 0;
    
    X_old = X_new;
    P_old = P_new;
    
    switch ii
        case 1 
            U = Bobs;
        otherwise
            U = conv2(X_old,P_old,'same');
    end
    divider = max(exp(-U/2),4.9407e-323);
    karmast = divider;
    etast = exp(-(Bobs-c).^2/(2*sigma2)).*U.*divider; 
    qstar = n_star(U,Bobs,c,sigma2); 

    q_ub = round(qstar + 5*sigma2.^.5); % m-n
    q_max = max(q_ub(:));

    for qs = 1:q_max
        switch qs
            case 700
                karmast = karmast.*U./qs;
                karmast700 = karmast.*divider;
                karmast = karmast700.*U./qs;
                keseist = karmast.*exp(-(Bobs-qs - c).^2/(2*sigma2));
                Sst = Sst.*divider+keseist*qs;
                etast = etast.*divider+keseist;
            otherwise
                karmast = karmast.*U./qs;
                keseist = karmast.*exp(-(Bobs-qs - c).^2/(2*sigma2));
                Sst = Sst+keseist*qs;
                etast = etast+keseist;
        end
    end
    
    EQst = Sst./etast;
    EQst(isnan(EQst)) = eps;
%% 
    H = psf2otf(P_new,size(Bobs));
    scale = real(ifftn(conj(H).*fftn(ones(size(H))))) + sqrt(eps);
    tmp = EQst./(ifftn(fftn(X_old).*H)+eps);
    X_new = X_old.*real(ifftn(conj(H).*fftn(tmp)))./scale;  
    HF_Reg = HessianFrobenious(X_old);
    X_new = X_new./(1+lambda_HF*HF_Reg);
    X_new = max(X_new,eps);
    
%%
    H = fftn(X_new);
    scale = otf2psf(conj(H).*fftn(ones(size(H))),size(P_old)) + sqrt(eps);
    tmp = EQst./(ifftn(H.*psf2otf(P_old,size(Bobs)))+eps);
    P_new = max(P_old.*(otf2psf(conj(H).*fftn(tmp),size(P_old)))./scale,0);
    P_new = P_new./(1+lambda_psf*P_old);
    P_new = AveragePSF(P_new);
    P_new = P_new/sum(P_new(:));
    if ~isfinite(P_new)
        disp('The psf is infinite!');
        break;
    end
    
%% 
    fun_all(ii,1) = norm(X_new - X_old);
    fun_all(ii,2) = norm(X_new);
    if ii>1
        fun_all(ii,3) = fun_all(ii,1)/fun_all(ii-1,2);
        if fun_all(ii,3) < tol
            break;
        end
    end
    
end
% fun_all((ii+1):end,:) = [];
X_out = X_new;
P_out = P_new;
end
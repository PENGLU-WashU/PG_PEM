function res = n_star(a,b,c,sigma2)
    res = sigma2/1.^2*lambertw(1/sigma2.*a.*exp(1/sigma2*(b-c)));
end
function HF = HessianFrobenious(Img)
    [Fx,Fy] = gradient(Img);
    [Fxx,~] = gradient(Fx);
    [Fxy,Fyy] = gradient(Fy);
    
    Hessian_norm = sqrt(Fxx.^2+2.*Fxy.^2+Fyy.^2);
    Fxx = Fxx./Hessian_norm;
    Fxy = Fxy./Hessian_norm;
    Fyy = Fyy./Hessian_norm;
    
    [Fxx_x,~] = gradient(Fxx);
    [Fxx_xx,~] = gradient(Fxx_x);
    
    [Fxy_x,~] = gradient(Fxy);
    [~,Fxy_xy] = gradient(Fxy_x);
    
    [~,Fyy_y] = gradient(Fyy);
    [~,Fyy_yy] = gradient(Fyy_y);
    
    HF = Fxx_xx+2*Fxy_xy+Fyy_yy;
end
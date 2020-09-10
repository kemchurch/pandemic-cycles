function [SI,p] = newton_F(SI,p,alpha,para,N,nu)
iter = 1;
while norm_nu_XbarN_vec(function_F(SI,p,alpha,para,N),nu,N)>1.8E-8 && iter<2000 &&  norm(function_F(SI,p,alpha,para,N))<1E4
    X = [SI;p] - function_DF(SI,p,alpha,para,N)\function_F(SI,p,alpha,para,N);
    SI = X(1:end-2);
    p = X(end-1:end);
    iter = iter+1;
end
end
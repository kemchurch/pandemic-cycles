function [SI,p,alpha,nF] = alpha_continuation(SI_0,p_0,alpha0,alpha1,steps,N,nu,para)
SI = zeros(8*(N+1),steps);
p = zeros(2,steps);
nF = zeros(1,steps);
alpha = linspace(alpha0,alpha1,steps);
[SI(:,1),p(:,1)] = newton_F(SI_0,p_0,alpha(1),para,N,nu);
nF(1) = norm_nu_XbarN_vec(function_F(SI(:,1),p(:,1),alpha(1),para,N),nu,N);
for n=2:steps
    [SI(:,n),p(:,n)] = newton_F(SI(:,n-1),p(:,n-1),alpha(n),para,N,nu);
    nF(n) = norm_nu_XbarN_vec(function_F(SI(:,n),p(:,n),alpha(n),para,N),nu,N);
end
end
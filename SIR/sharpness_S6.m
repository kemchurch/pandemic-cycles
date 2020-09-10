function pass = sharpness_S6(SI_0,SI_1,p0,p1,r0,para,N,nu,M,W,plotpar)
if nargin == 9
    plotpar = 0;
end
mu = intval(para(4));
gamma = intval(para(5));
w = intval(para(6));
SI_0 = intval(SI_0);
SI_1 = intval(SI_1);
p0 = intval(p0);
p1 = intval(p1);
nu = intval(nu);
SI_bar = SI_0 + infsup(0,1)*(SI_1-SI_0);
a_bar = intval(zeros(2*(N+1),4));
b_bar = intval(zeros(2*(N+1),4));
ab_bar = intval(zeros(2*(N+1),4));
p_bar = p0 + infsup(0,1)*(p1-p0);
eta = intval(zeros(2*(N+1),4));
Delta_eta = intval(zeros(1,4));
beta = [intval(para(2)),intval(para(2)),intval(para(3)),intval(para(3))];
z = [p_bar(1)-w, w, p_bar(2)-w-p_bar(1), w]/2;
abs_z = abs(z);
Delta_z = [r0/2,0,r0,0]/W(3);
for k=0:3
    a_bar(:,k+1) = intval([SI_bar(1+2*k*(N+1):(1+2*k)*(N+1));zeros(N+1,1)]);
    b_bar(:,k+1) = intval([SI_bar(1+(2*k+1)*(N+1):(2*k+2)*(N+1));zeros(N+1,1)]);
    ab_bar(:,k+1) = quadratic_cheb_intval(a_bar(:,k+1),b_bar(:,k+1));
    eta(:,k+1) = beta(k+1)*ab_bar(:,k+1) - (mu+gamma)*b_bar(:,k+1);
    Delta_eta(k+1) = abs_z(k+1)*(beta(k+1)*( norm_1_nu(a_bar(:,k+1),nu)/W(2) + norm_1_nu(b_bar(:,k+1),nu)/W(1) + r0/W(1)/W(2) ) + (mu+gamma)/W(2))*r0...
        + Delta_z(k+1)*( beta(k+1)*( norm_1_nu(ab_bar(:,k+1),nu) + r0*(norm_1_nu(a_bar(:,k+1),nu)/W(2) + norm_1_nu(b_bar(:,k+1),nu)/W(1)) + r0^2/W(1)/W(2)) + (mu+gamma)*(norm_1_nu(b_bar(:,k+1),nu) + r0/W(2)));
end
t = linspace(-1,1,M+1);
if plotpar==1
    tp = intval(zeros(M,1));
end
I_hat_prime = intval(zeros(M,4));
for m=1:M
    I_hat_prime(m,:) = z.*(eta(1,:) + 2*sum(eta(2:end,:).*repmat(ChebyshevT(1:2*N+1,infsup(t(m),t(m+1))).',[1,4]))) + infsup(-1,1)*Delta_eta;
    if plotpar==1
        tp(m)=infsup(t(m),t(m+1));
    end
end
check0 = min(I_hat_prime(:,1));
check1 = min(I_hat_prime(:,2));
check2 = min(I_hat_prime(:,3));
check3 = min(I_hat_prime(:,4));
if check0>0 & check1>0 & check2<0 & check3<0
    pass = 1;
else
    pass = 0;
end
if plotpar == 1
    plotintval([tp,I_hat_prime(:,1)].','r')
    hold on
    plotintval([tp,I_hat_prime(:,2)].','k')
    hold on
    plotintval([tp,I_hat_prime(:,3)].','b')
    hold on
    plotintval([tp,I_hat_prime(:,4)].','g')
    axis tight
end
end
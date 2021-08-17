function [Y0,Z0,Z1,Z2,radpol,radpol_fun,roots_double] = radpol_bounds(SI_0,SI_1,p_0,p_1,alpha_0,alpha_1,N,para,nu,W,rstar)
%para = (Lambda,beta0,beta1,mu,gamma,w,I_C,I_R)
nu = intval(nu);
rstar = intval(rstar);
alpha_0 = intval(alpha_0);
alpha_1 = intval(alpha_1);
SI_0 = intval(SI_0);
SI_1 = intval(SI_1);
p_0 = intval(p_0);
p_1 = intval(p_1);
A_dagger = function_DF_int(SI_0,p_0,alpha_0,para,N);
beta_k = abs([intval(para(2)),intval(para(2)),intval(para(3)),intval(para(3))]);
w = intval(para(6));
mu = intval(para(4));
gamma = intval(para(5));
abs_z_k = abs([p_0(1)-w, w, p_0(2)-w-p_0(1), w])/2;
delta_SI = SI_1 - SI_0;
delta_p = p_1 - p_0;
delta_alpha = alpha_1 - alpha_0;
Delta = [delta_SI;delta_p;delta_alpha];

% Bar calculations
Abar = inv(A_dagger);
SI_bars = SI_0 + infsup(0,1)*(SI_1-SI_0);
p_bars = p_0 + infsup(0,1)*(p_1-p_0);
alpha_bars = alpha_0 + infsup(0,1)*(alpha_1-alpha_0);
A_dagger_bars_ext = [function_DF_int(SI_bars,p_bars,alpha_bars,para,N),...
    function_DF_alpha_int(SI_bars,alpha_bars,para,N)];

% Norm bar calculations
nAFbar = norm_nu_XbarN_vec(Abar*A_dagger_bars_ext*Delta,nu,N,W); % INCOMPLETE.
nAF0 = norm_nu_XbarN_vec(Abar*function_F_int(SI_0,p_0,alpha_0,para,N),nu,N,W);

% Precomputation: Y0
TPsi = intval(zeros(2*2*(N+1),4));
nn = 1:1:2*N+1;
wY0 = intval(zeros(2*N+2));
wY0(2:end,2:end) = diag(nu.^(nn)/2./nn);
SI_bars_k = intval(zeros(2*(N+1),4));
for k=0:3
    SI_bars_k(:,k+1) = SI_bars(1+2*k*(N+1):2*(k+1)*(N+1));
    TPsi(:,k+1) = abs(function_TPsi2_int(SI_bars_k(1:N+1,k+1),SI_bars_k(N+2:end,k+1),p_bars,para,k,N));
    TPsi(1:N,k+1) = 0;
    TPsi(2*(N+1)+1:2*(N+1)+N,k+1) = 0;
    TPsi(:,k+1) = repmat(wY0,[2,2])*TPsi(:,k+1);
end

% Bound Y0 -- note, we take the max of W(1) and W(2) for the weight on the
% TPsi terms because the latter are already very small and have almost
% negligible impact on the Y0 bound.
Y0 = nAFbar + nAF0 + max(W(1),W(2))*max(sum(TPsi,1));

% Bound Z0
I_XbarN = intval(eye(8*(N+1)+2,8*(N+1)+2));
Z0 = norm_nu_XbarN_mat_int(I_XbarN-Abar*A_dagger,nu,N,W);

% Precomputation: Z1
norm_gz = [1/intval(2),0,1,0];
Z11k = intval(zeros(2*(N+1),4));
Z12k = intval(zeros(2,4));
Psi = intval(zeros(2*(N+1),8));
na_k = intval(zeros(4,1));
nb_k = intval(zeros(4,1));
lambda = [mu;mu+gamma];
for k=0:3
    Z11k(:,k+1) = beta_k(k+1)*abs_z_k(k+1)*repmat(convcheb_tailbound(SI_0(1+2*k*(N+1):(2*k+1)*(N+1)),nu)/W(2)...
        + convcheb_tailbound(SI_0(1+(2*k+1)*(N+1):(2*k+2)*(N+1)),nu)/W(1), [2,1]);
    Psi(:,1+2*k:2*(k+1)) = function_TPsi2_int(SI_0(1+2*k*(N+1):(2*k+1)*(N+1)),SI_0(1+(2*k+1)*(N+1):(2*k+2)*(N+1)),...
        p_0,para,k,N,'no_T');
    na_k(k+1)=norm_1_nu(SI_bars_k(1:N+1,k+1),nu);
    nb_k(k+1)=norm_1_nu(SI_bars_k(N+2:end,k+1),nu);
    Z12k(:,k+1)=W(1:2).*( beta_k(k+1)*(na_k(k+1)/W(2) + nb_k(k+1)/W(1))*[1;1] + beta_k(k+1)*(na_k(k+1)/W(2) + nb_k(k+1)/W(1))*[1;1] + lambda./W(1:2)+ norm_gz(k+1)*norm_1_nu(Psi(:,1+2*k:2*(k+1)),nu).' );
end
bold_h = [reshape(Z11k,[8*(N+1),1]);0;0]*infsup(-1,1);
bold_W = diag([repmat([W(1)*ones(1,N+1),W(2)*ones(1,N+1)],[1,4]),W(3),W(3)]);
bold_W_inv = inv(bold_W);
Hbar = intval(zeros(8*(N+1)+2,1));
rbar = intval(zeros(8*(N+1)+2,1));
Gbar = intval(zeros(8*(N+1)+2,1));
Hbar(1:N+1:8*(N+1))=infsup(-1,1);
rbar(N+1:N+1:8*(N+1))=infsup(-1,1);
Gbar([end-1,end])=infsup(-1,1);
nT = (2*nu+1/nu);

% Bound Z1
Z1 = nT*norm_nu_XbarN_vec(abs(Abar)*bold_h,nu,N,W)...
    + nT/2/(N+1)*max(max(Z12k))...
    + 1/nu^(N+1)*(norm_nu_XbarN_vec(Abar*bold_W_inv*Hbar,nu,N,W) + norm_nu_XbarN_vec(Abar*bold_W_inv*Gbar,nu,N,W) + (mu+gamma)*max(abs_z_k)/2*norm_nu_XbarN_vec(Abar*bold_W_inv*rbar,nu,N,W));

% Precomputation: Z2
abs_gz = intval([1/2,0;0,0;1/2,1/2;0,0]);
xi_hat_1 = intval(zeros(4,2,3));
xi_hat_2 = intval(zeros(4,2,3));
n_diff_ak = intval(zeros(4,1));
n_diff_bk = intval(zeros(4,1));
n_diffa_diffb = intval(zeros(4,1));
n_a_diffb = intval(zeros(4,1));
n_b_diffa = intval(zeros(4,1));
for k=0:3
     n_diffa_diffb(k+1) = norm_1_nu(quadratic_cheb_intval(delta_SI(1+2*k*(N+1):(2*k+1)*(N+1)),delta_SI(1+(2*k+1)*(N+1):(2*k+2)*(N+1))), nu);
     n_a_diffb(k+1) = norm_1_nu(quadratic_cheb_intval(SI_0(1+2*k*(N+1):(1+2*k)*(N+1)),delta_SI(1+(2*k+1)*(N+1):(2*k+2)*(N+1))), nu);
     n_b_diffa(k+1) = norm_1_nu(quadratic_cheb_intval(SI_0(1+(2*k+1)*(N+1):(2*k+2)*(N+1)),delta_SI(1+2*k*(N+1):(1+2*k)*(N+1))), nu);
     n_diff_ak(k+1) = norm_1_nu(delta_SI(1+2*k*(N+1):(2*k+1)*(N+1)),nu);
     n_diff_bk(k+1) = norm_1_nu(delta_SI(1+(2*k+1)*(N+1):(2*k+2)*(N+1)),nu);
    for j=1:2
        % quadratic terms
        xi_hat_1(k+1,j,1) = 1/W(3)*norm_gz(k+1)*beta_k(k+1)*2/W(1)/W(2);
        xi_hat_2(k+1,j,1) = beta_k(k+1)/W(1)/W(2);
        % linear terms
        xi_hat_1(k+1,j,2) = abs_z_k(k+1)*beta_k(k+1)*(2/W(2)/W(1)) + abs_gz(k+1,:)*delta_p*(beta_k(k+1)*2/W(1)/W(2)) + norm_gz(k+1)/W(3)*(beta_k(k+1)*(1/W(2)*(na_k(k+1)+n_diff_ak(k+1)) + 1/W(1)*(nb_k(k+1)+n_diff_bk(k+1)) ) + lambda(j)/W(j) );
        xi_hat_2(k+1,j,2) = beta_k(k+1)*(1/W(1)*(nb_k(k+1)+n_diff_bk(k+1)) + 1/W(2)*(na_k(k+1)+n_diff_ak(k+1))) + lambda(j)/W(j);
        % constant terms
        xi_hat_1(k+1,j,3) = abs_z_k(k+1)*beta_k(k+1)*(1/W(2)*n_diff_ak(k+1) + 1/W(1)*n_diff_bk(k+1)) + abs_gz(k+1,:)*delta_p*(beta_k(k+1)*(1/W(2)*(na_k(k+1)+n_diff_ak(k+1)) + 1/W(1)*(nb_k(k+1) + n_diff_bk(k+1)) ) + lambda(j)/W(j) ) ;
        xi_hat_2(k+1,j,3) = beta_k(k+1)*(n_diffa_diffb(k+1)+n_a_diffb(k+1)+n_b_diffa(k+1));
    end
end
xi_hat_1_rstar = intval(zeros(4,2));
xi_hat_2_rstar = intval(zeros(4,2));
for k=0:3
    for j=1:2
        xi_hat_1_rstar(k+1,j) = reshape(xi_hat_1(k+1,j,:),[1,3])*[rstar^2;rstar;1];
        xi_hat_2_rstar(k+1,j) = reshape(xi_hat_2(k+1,j,:),[1,3])*[rstar^2;rstar;1];
    end
end
% The A_bold terms
A11 = Abar(1:end-2,1:end-2);
A12 = Abar(1:end-2,end-1:end);
A21 = Abar(end-1:end,1:end-2);
nA11 = intval(zeros(4,2,4,2));
nA12 = intval(zeros(4,2));
nA21 = intval(zeros(2,4,2));
nA22 = sum(abs(Abar(end-1:end,end-1:end)),2);
% Build nA11
for k=0:3
    for j=1:2
        for m=0:3
            for ell=1:2
                i = 2*k+j-1;
                q = 2*m+ell-1;
                nA11(k+1,j,m+1,ell) = norm_nu_1_op(A11(1+i*(N+1):(i+1)*(N+1), 1+q*(N+1):(q+1)*(N+1)),nu);
            end
        end
    end
end
% Build nA12
for m=0:3
    for ell=1:2
        q = 2*m+ell-1;
        nA12(m+1,ell) = norm_1_nu(A12(1+q*(N+1):(q+1)*(N+1),end-1).',nu) + norm_1_nu(A12(1+q*(N+1):(q+1)*(N+1),end).',nu);
    end
end
% Build nA21
for j=1:2
    for m=0:3
        for ell=1:2
            q = 2*m+ell-1;
            nA21(j,m+1,ell) = norm_1_omega_dual(A21(end-1+(j-1),1+q*(N+1):(q+1)*(N+1)),nu);
        end
    end
end
% Build mu_hat and mu_hat_inf
mu_hat = intval(zeros(4,2));
mu_hat_inf = intval(zeros(1,2));
for j=1:2
    mu_hat_inf(j) = sum(sum(reshape(nA21(j,:,:),[4,2]).*(xi_hat_1_rstar+repmat(norm_gz.',[1,2]).*xi_hat_2_rstar),1)) + nA22(j)/W(2)*delta_alpha + delta_alpha/W(2)/(2*(N+1));
    for k=0:3
        mu_hat(k+1,j) = sum(sum(reshape(nA11(k+1,j,:,:),[4,2]).*(xi_hat_1_rstar+repmat(norm_gz.',[1,2]).*xi_hat_2_rstar),1)) + nA12(k+1,j)/W(2)*delta_alpha + (xi_hat_1_rstar(k+1,j)+norm_gz(k+1)*xi_hat_2_rstar(k+1,j)/W(3))/(2*(N+1));
    end
end
% Z2 bound
Z2 = max([reshape(repmat(W(1:2)',[4,1]).*mu_hat, [1,8]) , W(3)*mu_hat_inf]);

% Radii polynomial
radpol = [Z2+Z1+Z0-1,Y0];
roots_double = sort(roots(sup(radpol)));
radpol_fun=@(x)radpol(1)*x + radpol(2);
end
function F = cheb_ODE(S,I,Splus,Iplus,p,para,idx,N,dflag)
% Notes:
% para = (Lambda,beta0,beta1,mu,gamma,w,I_C,I_R)
if nargin==8
    dflag=1;
end
% Note: the dflag parameter is used so that this function can be recycled
% in the partial derivative calculation

% Parameters
Lambda = para(1);
if idx==0 || idx==1
    beta = para(2);
elseif idx==2 || idx==3
    beta = para(3);
end
mu = para(4);
gamma = para(5);
w = para(6);
if idx==0
    z = (p(1)-w)/2;
elseif idx==1
    z = w/2;
elseif idx==2
    z = (p(2)-w-p(1))/2;
elseif idx==3
    z = w/2;
end

% Preallocation
n = ones(1,N+1);
T = -diag(n,-1) + diag(n,1);
T(1,:) = 0;
SI_conv = quadratic_cheb([S;0],[I;0]);
F = dflag*repmat(2*(0:N)',[2,1]).*[S;I];
Id = eye(N+2);
Sdot = z*(Lambda*Id(:,1) - beta*SI_conv - mu*[S;0]);
Idot = z*(beta*SI_conv - (mu+gamma)*[I;0]);
TSdot = T*Sdot;
TIdot = T*Idot;
alt = 1-2*rem(1:N,2);
S_zero = S(1)+2*sum(S(2:end)) - Splus(1) - 2*sum(alt'.*Splus(2:end));
I_zero = I(1)+2*sum(I(2:end)) - Iplus(1) - 2*sum(alt'.*Iplus(2:end));

% Equations
F(1) = dflag*S_zero;
F(N+2) = dflag*I_zero;
F(2:N+1) = F(2:N+1) + TSdot(2:N+1);
F(N+3:end) = F(N+3:end) + TIdot(2:N+1);

end
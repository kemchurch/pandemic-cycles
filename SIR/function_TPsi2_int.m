function Psi = function_TPsi2_int(S,I,p,para,idx,N,Tconv)
% Notes:
% para = (Lambda,beta0,beta1,mu,r,w,Iminus,Iplus)

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
    c = (p(1)-w)/2;
elseif idx==1
    c = w/2;
elseif idx==2
    c = (p(2)-w-p(1))/2;
elseif idx==3
    c = w/2;
end

if nargin==6
    Tconv=1;
end

% Preallocation
n = ones(1,2*(N+1)-1);
T = -diag(n,-1) + diag(n,1);
T(1,:) = 0;
S_ext = [S;intval(zeros(N+1,1))];
I_ext = [I;intval(zeros(N+1,1))];
SI_conv = quadratic_cheb_intval(S_ext,I_ext);
Id = intval(eye(2*(N+1)));
Sdot = c*(Lambda*Id(:,1) - beta*SI_conv - mu*S_ext);
Idot = c*(beta*SI_conv - (mu+gamma)*I_ext);
if Tconv==1
    TSdot = T*Sdot;
    TIdot = T*Idot;
end

% Equations
if Tconv==1
    Psi = [0;TSdot(2:end);0;TIdot(2:end)];
else
    Psi = [Sdot,Idot];
end

end
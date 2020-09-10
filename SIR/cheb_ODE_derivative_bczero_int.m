function [DFS,DFI] = cheb_ODE_derivative_bczero_int(S,I,p,para,idx,N)
% Parameters
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
n = (1:N+2)';
J = (1:N+1).*ones(N+2,N+1);
Id = eye(N+2);
nminusJ = abs(n-J)+1;
nplusJ = abs(n+J)-1;
DFS = intval([diag(2*(0:N)'),zeros(N+1,N+1)]);
DFI = intval([zeros(N+1,N+1),diag(2*(0:N)')]);
I_ext = intval([I;zeros(N+2,1)]);
S_ext = intval([S;zeros(N+2,1)]);
DFS_nL = intval(zeros(N+2,2*(N+1)));
DFI_nL = intval(zeros(N+2,2*(N+1)));

% Compute convolution parts
DFS_nL(:,1) = z*( -mu*Id(:,1) - beta*I_ext(nminusJ(:,1)) );
DFS_nL(:,N+2) = z*( -beta*S_ext(nminusJ(:,1)) );
DFS_nL(:,2:N+1) = z*( -mu*Id(:,2:N+1) - beta*( I_ext(nminusJ(:,2:end)) + I_ext(nplusJ(:,2:end)) ) );
DFS_nL(:,N+3:end) = z*( -beta*( S_ext(nminusJ(:,2:end)) + S_ext(nplusJ(:,2:end)) ) );
DFI_nL(:,1) = z*( beta*I_ext(nminusJ(:,1)) );
DFI_nL(:,N+2) = z*( beta*S_ext(nminusJ(:,1)) - (mu+gamma)*Id(:,1) );
DFI_nL(:,2:N+1) = z*( beta*( I_ext(nminusJ(:,2:end)) + I_ext(nplusJ(:,2:end)) ) );
DFI_nL(:,N+3:end) = z*( beta*( S_ext(nminusJ(:,2:end)) + S_ext(nplusJ(:,2:end)) ) - (mu+gamma)*Id(:,2:N+1) );
TDFS_nL = T*DFS_nL;
TDFI_nL = T*DFI_nL;
% Outputs
DFS = DFS + TDFS_nL(1:N+1,:);
DFI = DFI + TDFI_nL(1:N+1,:);

end
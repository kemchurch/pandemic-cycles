function nrm = norm_nu_XbarN_vec(x,nu,N,W)
%Note: Norm in X
if nargin == 3
    W = [1;1;1];
end
xr = reshape(x(1:8*(N+1)),[N+1,8])*diag(repmat(W(1:2),[4,1]));
cp = x(end-1:end).*W(3);
norm_xr_nu = norm_1_nu(xr,nu);
norm_cp = norm(cp,inf);
nrm = max([norm_xr_nu,norm_cp]);
end
function F = function_F(SI,p,alpha,para,N)
% Notes:
% para = (Lambda,beta0,beta1,mu,gamma,w,I_C,I_R)
% s = convexity parameter
[S,I]=convert_1to2(SI,N);
S = reshape(S,[N+1,4]);
I = reshape(I,[N+1,4]);
F = [cheb_ODE(S(:,1),I(:,1),S(:,2),I(:,2),p,para,0,N);
    cheb_ODE(S(:,2),I(:,2),S(:,3),I(:,3),p,para,1,N);
    cheb_ODE(S(:,3),I(:,3),S(:,4),I(:,4),p,para,2,N);
    cheb_ODE(S(:,4),I(:,4),S(:,1),I(:,1),p,para,3,N);
    cheb_threshold(I(:,2),para,alpha,1,N);
    cheb_threshold(I(:,4),para,alpha,3,N)];
end
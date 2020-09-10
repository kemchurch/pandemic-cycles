function F = function_F_int(SI,p,alpha,para,N)
% Notes:
% para = (Lambda,beta0,beta1,mu,r,w,Iminus,Iplus)
% s = convexity parameter
[S,I]=convert_1to2_int(SI,N);
S = reshape(S,[N+1,4]);
I = reshape(I,[N+1,4]);
F = [cheb_ODE_int(S(:,1),I(:,1),S(:,2),I(:,2),p,para,0,N);
    cheb_ODE_int(S(:,2),I(:,2),S(:,3),I(:,3),p,para,1,N);
    cheb_ODE_int(S(:,3),I(:,3),S(:,4),I(:,4),p,para,2,N);
    cheb_ODE_int(S(:,4),I(:,4),S(:,1),I(:,1),p,para,3,N);
    cheb_threshold_int(I(:,2),para,alpha,1,N);
    cheb_threshold_int(I(:,4),para,alpha,3,N)];
end
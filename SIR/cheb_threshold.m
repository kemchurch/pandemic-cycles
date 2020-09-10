function thr = cheb_threshold(I,para,alpha,idx,N)
% Notes:
% para = (Lambda,beta0,beta1,mu,r,w,I_C,I_R)
% s = convexity parameter

% Parameters
if idx==1
    Icrit = para(7);
elseif idx==3
    Icrit = para(8);
end

% Preallocation
alt = zeros(1,N+1);
alt(3:2:N+1) = 2;
alt(3:N+1) = alt(3:N+1)./((2:N).^2-1);

% Output
thr = Icrit - alpha*(I(1)-sum(alt'.*I)) - (1-alpha)*(I(1)+2*sum(I(2:N+1)));
end
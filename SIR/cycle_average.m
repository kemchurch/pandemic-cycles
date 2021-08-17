function [avg,radius] = cycle_average(SI,p,para,N,nu,r0,W)
[S,I]=convert_1to2(SI,N);
I = reshape(intval(I),[N+1,4]);
S = reshape(intval(S),[N+1,4]);
p = intval(p);
r0 = intval(r0);
alt = intval([1,zeros(1,N)]);
alt(3:2:end)=2;
alt(3:end) = alt(3:end)./(intval(2:N).^2-1);
w = intval(para(6));
z = [p(1)-w, w, p(2)-w-p(1), w]/2;
for n=1:N+1
    I(n,:) = I(n,:).*z;
    S(n,:) = S(n,:).*z;
end
avg_intval = (2/p(2))*sum(alt*I);
avg = mid(avg_intval);
Delta_z = [r0/2,0,r0,0]/W(3);
err = intval(zeros(4,1));
for k=0:3
    err(k+1) = (p(2)*Delta_z(k+1) + r0/W(3)*z(k+1))*(r0/W(2) + norm_1_nu(I(:,k+1),nu))/(p(2)*(p(2)-r0/W(3)));
end
radius = sup(r0/W(2) + sum(err)) + sup(abs(avg - avg_intval));
end
function [min,max,radius] = cycle_minmax(SI,p,N,r0,W)
[S,I]=convert_1to2(SI,N);
I = reshape(intval(I),[N+1,4]);
S = reshape(intval(S),[N+1,4]);
p = intval(p);
r0 = intval(r0);
radius = sup(r0/W(2));
min = inf(I(1,4) + 2*sum(I(2:end,4)));
max = sup(I(1,2) + 2*sum(I(2:end,2)));
end
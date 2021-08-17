function pass = sharpness_S5(SI_0,SI_1,r0,N,para,W)
I_C = intval(para(7));
I_R = intval(para(8));
b0 = SI_0(N+2:2*(N+1)) + infsup(0,1)*(SI_1(N+2:2*(N+1)) - SI_0(N+2:2*(N+1)));
b2 = SI_0(1+5*(N+1):6*(N+1)) + infsup(0,1)*(SI_1(1+5*(N+1):6*(N+1)) - SI_0(1+5*(N+1):6*(N+1)));
check0 = b0(1) + 2*sum(b0(2:end));
check2 = b2(1) + 2*sum(b2(2:end));
if I_R+r0/W(2) < check0 & check2 < I_C-r0/W(2)
    pass = 1;
else
    pass = 0;
end
end
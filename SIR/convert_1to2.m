function [S,I] = convert_1to2(SI,N)
S = zeros(4*(N+1),1);
I = zeros(4*(N+1),1);
for k=0:3
    S(1+k*(N+1):(k+1)*(N+1)) = SI(1+k*2*(N+1):k*2*(N+1)+N+1);
    I(1+k*(N+1):(k+1)*(N+1)) = SI(1+k*2*(N+1)+N+1:(2*k+2)*(N+1));
end
end


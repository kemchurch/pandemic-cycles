function DF = function_DF_int(SI,p,alpha,para,N)
[S,I]=convert_1to2_int(SI,N);
S = reshape(S,[N+1,4]);
I = reshape(I,[N+1,4]);
DF = intval(zeros(2+4*2*(N+1),2+4*2*(N+1)));
alt = 1-2*rem(1:N,2);

% Derivatives of BVPs wrt (S,I)_idx
for idx=0:3
    % n/=/0
    [DFS0,DFI0] = cheb_ODE_derivative_bczero_int(S(:,idx+1),I(:,idx+1),p,para,idx,N);
    DF(1+idx*2*(N+1):idx*2*(N+1)+N+1,1+idx*2*(N+1):(idx+1)*2*(N+1)) = DFS0; 
    DF(N+2+idx*2*(N+1):2*(N+1)+idx*2*(N+1),1+idx*2*(N+1):2*(N+1)+idx*2*(N+1)) = DFI0; 
    % n=0 --
    midx = mod(idx+1,4);
    DF(1+idx*2*(N+1),1+idx*2*(N+1):idx*2*(N+1)+N+1) = [1,2*ones(1,N)]; 
    DF(1+idx*2*(N+1),1+midx*2*(N+1):midx*2*(N+1)+N+1) = [-1,-2*alt]; 
    DF(N+2+idx*2*(N+1),N+2+idx*2*(N+1):2*(N+1)+idx*2*(N+1)) = [1,2*ones(1,N)]; 
    DF(N+2+idx*2*(N+1),N+2+midx*2*(N+1):2*(N+1)+midx*2*(N+1)) = [-1,-2*alt]; 
end

% Derivatives wrt p
DF(:,end-1:end) = cheb_ODE_derivative_bczero_para_int(S,I,para,N);

% Derivatives of thresholds wrt (S,I)_idx
D = cheb_threshold_derivative(alpha,N);
DF(end-1,1+3*(N+1):2*(N+1)+2*(N+1)) = D;
DF(end,end-(N+2):end-2) = D;
end
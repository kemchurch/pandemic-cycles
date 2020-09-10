function nrm = norm_nu_XbarN_mat_int(A,nu,N,W)
% Note: linear operator X->X
nu = intval(nu);
A11 = A(1:end-2,1:end-2);
A12 = A(1:end-2,end-1:end);
A21 = A(end-1:end,1:end-2);
A22 = A(end-1:end,end-1:end);
omega = [1,2*nu.^(1:1:N)];
normsA11 = intval(zeros(8,1));
normsA12 = intval(zeros(8,2));
normsA21 = intval(zeros(2,8));
W12 = W(1:2);
W12rep = repmat(W12,[4,1]);

for i=0:7
    normsA12(i+1,:) = norm_1_nu(A12(1+i*(N+1):(i+1)*(N+1),:),nu)*W(1+mod(i,2));
    normsA21(:,i+1) = sum(abs(A21(:,1+i*(N+1):(i+1)*(N+1)))./omega./W(1:2),2);
    for q=0:7
        normsA11(i+1) = normsA11(i+1) + norm_nu_1_op(A11(1+i*(N+1):(i+1)*(N+1), 1+q*(N+1):(q+1)*(N+1)),nu,W12(1+mod(q+1,2)));
    end
end
nA11 = max(normsA11.*W12rep);
nA22 = norm(A22,inf);
nA12 = (max(normsA12(:,1)) + max(normsA12(:,2)))/W(3);
nA21 = max(max(normsA21(1,:)), max(normsA21(2,:)) )*W(3);
nrm = nA11 + nA12 + nA21 + nA22;
end
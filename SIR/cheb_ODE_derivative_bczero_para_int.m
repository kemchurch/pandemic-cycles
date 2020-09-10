function Dpara = cheb_ODE_derivative_bczero_para_int(S,I,para,N)
% para = (Lambda,beta0,beta1,mu,r,w,Iminus,Iplus)
% Note: a hack for this is to simply set the w parameter to zero and
% set the p parameters to either 1 or zero, depending on which partial
% derivative is needed. Then, copy the F map with dflag=0 to get rid of the 
% diagonal-dominant terms and the zero rows that do not depend on the
% parameter. Then just kill off the rows for which we know the derivative 
% is zero anyway.
S = reshape(S,[N+1,4]);
I = reshape(I,[N+1,4]);
% Adjust the w parameter to zero:
para(6)=0;

Fp0 = [cheb_ODE_int(S(:,1),I(:,1),S(:,2),I(:,2),[1,0],para,0,N,0);
    intval(zeros(2*(N+1),1));
    cheb_ODE_int(S(:,3),I(:,3),S(:,4),I(:,4),[1,0],para,2,N,0);
    intval(zeros(2*(N+1),1));
    intval(0);
    intval(0)];
Fp1 = [intval(zeros(2*(N+1),1));
    intval(zeros(2*(N+1),1));
    cheb_ODE_int(S(:,3),I(:,3),S(:,4),I(:,4),[0,1],para,2,N,0);
    zeros(2*(N+1),1);
    intval(0);
    intval(0)];
Dpara = [Fp0,Fp1];

end
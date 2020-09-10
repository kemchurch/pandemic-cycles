function [pass,radii] = verify_branch(SI,times,alpha,para,N,M,nu,parallel,W,rstar)
% Running in parallel: parallel = [x,y] means load will be split into y
% computations, with x\in{1,2,..,y} denoting the current thread.
if isempty(parallel)
    parallel = [1,1];
end
branches = length(alpha);
bsplit = floor(branches/parallel(2));
n1 = 1+(parallel(1)-1)*bsplit;
if parallel(1)==parallel(2)
    n2 = branches;
else
    n2 = parallel(1)*bsplit+1;
end
pass = 1;
radii = intval(zeros(n2-n1,1));
for n=n1:n2-1
   [~,~,~,~,~,radfun,root] = radpol_bounds(SI(:,n),SI(:,n+1),times(:,n),times(:,n+1),alpha(n),alpha(n+1),N,para,nu,W,rstar);
   if root(1)<=0
       disp(['Radii polynomial failed at n=',num2str(n)])
       pass = 0;
       break
   end
   if root(1)>=rstar
       disp(['Failure with rstar at n=',num2str(n)])
       pass = 0;
       break
   end
   iter = 1;
   r0 = intval(root(1));
   while not(radfun(r0)<0) & r0<rstar & iter < 100
       r0 = r0 + eps;
       iter = iter+1;
   end
   if r0>=rstar | radfun(r0)>=0
       disp(['Something is wrong with the roots at n=',num2str(n)])
       pass = 0;
       break
   end
   S4 = sharpness_S4(times(:,n),times(:,n+1),para);
   S5 = sharpness_S5(SI(:,n),SI(:,n+1),r0,N,para,W);
   S6 = sharpness_S6(SI(:,n),SI(:,n+1),times(:,n),times(:,n+1),r0,para,N,nu,M,W,0);
   if min([S4,S5,S6])==0
       disp(['A sharpness check failed at n=',num2str(n)])
       pass = 0;
       break
   else
       disp(['Success at n=',num2str(n)])
       radii(n+1-n1) = r0;
   end
end
if pass==1
    disp(['Success: process thread ',num2str(parallel(1)),' of ',num2str(parallel(2)),' proven successfully.']);
else
    disp('Failure.')
end
end
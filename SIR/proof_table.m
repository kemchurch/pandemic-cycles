% This script proves and computes the table from Theorem 12 and the 
% overshooting duration.
clear all
load('candidate_cheb_N16.mat')
TABLE(1,:) = ["alpha",0,"1"];
% alpha = 0 computations
[~,~,~,~,~,radpol_fun0,r0_alpha0] = radpol_bounds(SI(:,1),SI(:,1),times(:,1),times(:,1),alpha(1),alpha(1),N,para,nu,W,rstar);
if radpol_fun0(r0_alpha0)>=0 | r0_alpha0>rstar
    disp('Error.') %Returns an error if the radii polynomial fails.
    return
end
[I_avg_alpha0,I_avg_alpha0_err] = cycle_average(SI(:,1),times(:,1),para,N,nu,r0_alpha0,W);
c_alpha0 = times(1,1);
p_alpha0 = times(2,1);
c_err_alpha0 = r0_alpha0/W(3);
p_err_alpha0 = r0_alpha0/W(3);
Iplus_alpha0 = 1600;
Iplus_alpha0_err = 0;

% alpha = 1 computations
[~,~,~,~,~,radpol_fun1,r0_alpha1] = radpol_bounds(SI(:,end),SI(:,end),times(:,end),times(:,end),alpha(end),alpha(end),N,para,nu,W,rstar);
if radpol_fun1(r0_alpha1)>=0 | r0_alpha1>rstar
    disp('Error.') %Returns an error if the radii polynomial fails.
    return
end
[I_avg_alpha1,I_avg_alpha1_err] = cycle_average(SI(:,end),times(:,end),para,N,nu,r0_alpha1,W);
[~,Iplus_alpha1,Iplus_alpha1_err] = cycle_minmax(SI(:,end),times(:,end),N,r0_alpha1,W);
c_alpha1 = times(1,end);
p_alpha1 = times(2,end);
c_err_alpha1 = r0_alpha1/W(3);
p_err_alpha1 = r0_alpha1/W(3);

% Overshoot
[t0,t1,check] = overshoot(SI(:,end),times(:,end),N,para,r0_alpha1,W,1E-14);
if check == 0
    disp('Error') %Returns an error if the bisection/Newton procedure fails to find good t0 and t1.
    return
end
overshoot = mid(t1-t0);
overshoot_error = sup(abs(t1-t0-mid(t1-t0)));

% Build the table
TABLE(2,:) = ["c",num2str(c_alpha0,20),num2str(c_alpha1,20)];
TABLE(3,:) = ["c_err",num2str(c_err_alpha0,20),num2str(c_err_alpha1,20)];
TABLE(4,:) = ["p",num2str(p_alpha0,20),num2str(p_alpha1,20)];
TABLE(5,:) = ["p_err",num2str(p_err_alpha0,20),num2str(p_err_alpha1,20)];
TABLE(6,:) = ["I^plus",num2str(Iplus_alpha0,20),num2str(Iplus_alpha1,20)];
TABLE(7,:) = ["I^plus_err",num2str(Iplus_alpha0_err,20),num2str(Iplus_alpha1_err,20)];
TABLE(8,:) = ["[I]",num2str(I_avg_alpha0,20),num2str(I_avg_alpha1,20)];
TABLE(9,:) = ["[I]_err",num2str(I_avg_alpha0_err,20),num2str(I_avg_alpha1_err,20)];
TABLE(10,:) = ["r_0",num2str(r0_alpha0,20),num2str(r0_alpha1,20)];

% Display the output
disp(TABLE);
disp(['Overshoot duration: ',num2str(overshoot,20),'.']);
disp(['Overshoot duration error bound: ',num2str(overshoot_error,20),'.']);

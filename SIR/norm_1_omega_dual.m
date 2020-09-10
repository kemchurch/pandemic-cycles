function Norm = norm_1_omega_dual(a,nu)
N = length(a(:,1))-1;
omega = nu.^((0:N)');
omega(2:end) = 2*omega(2:end);
Norm = max(abs(a)./omega);
end
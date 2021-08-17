function D = cheb_threshold_derivative(alpha,N)
% Note: This is the derivative of the threshold with respect to the only
% variable, so it will need to be dropped into the appropriate place later.

% Preallocation
alt = zeros(1,N+1);
alt(3:2:N+1) = 2;
alt(3:N+1) = alt(3:N+1)./((2:N).^2-1);

% Output
D = - alpha*[1,-alt(2:end)] - (1-alpha)*[1,2*ones(1,N)];
end
function y = ChebyshevTprime(n,x)
y = n.*sin(n.*acos(x))/sqrt(1-x^2);
end
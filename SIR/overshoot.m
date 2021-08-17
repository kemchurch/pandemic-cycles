function [t0,t1,check] = overshoot(SI,p,N,para,r0,W,tol)
I_C = para(7);
[~,I]=convert_1to2(SI,N);
I = reshape(I,[N+1,4]);
t0bar=[0,0];
t1bar=[0,0];
% Get t0bar with double arithmetic
for k=0:1
    f = @(x)(-1)^k*(r0+tol)/W(2) - I_C +  I(1,2) + 2*ChebyshevT(1:N,x)*I(2:end,2);
    f_diff = @(x) 2*ChebyshevTprime(1:N,x)*I(2:end,2);
    % Crude bisection to get an initial guess for t0bar
    a=-1; b=1;
    fa = f(a);
    for n=1:300
        t = (a+b)/2;
        ft = f(t);
        if fa*ft>=0
            a = t;
        else
            b = t;
            fa = ft;
        end
    end
    % Refinement with Newton
    iter = 1;
    while abs(f(t))>1E-16 & iter<10000
        t = t - f(t)/f_diff(t);
        iter = iter+1;
    end
    t0bar(k+1)=real(t);
end
% Get t1bar with double arithmetic
for k=0:1
    f = @(x) -(-1)^k*(r0+tol)/W(2) - I_C +  I(1,3) + 2*ChebyshevT(1:N,x)*I(2:end,3);
    f_diff = @(x) 2*ChebyshevTprime(1:N,x)*I(2:end,3);
    % Crude bisection to get an initial guess for t1bar
    a=-1; b=1;
    fa = f(a);
    for n=1:300
        t = (a+b)/2;
        ft = f(t);
        if fa*ft>=0
            a = t;
        else
            b = t;
            fa = ft;
        end
    end
    % Refinement with Newton
    iter = 1;
    while abs(f(t))>1E-16 & iter<10000
        t = t - f(t)/f_diff(t);
        iter = iter+1;
    end
    t1bar(k+1)=real(t);
end
% Check with interval arithmetic that the computed t0bar and t1bar satisfy
% the inequalities the lemma.
r0 = intval(r0);
I = intval(I);
I_C = intval(I_C);
t0bar = intval(t0bar);
t1bar = intval(t1bar);
certificate = [r0/W(2) - I_C +  I(1,2) + 2*ChebyshevT(1:N,t0bar(1))*I(2:end,2);
    - (-r0/W(2) - I_C +  I(1,2) + 2*ChebyshevT(1:N,t0bar(2))*I(2:end,2));
    - ( -r0/W(2) - I_C +  I(1,3) + 2*ChebyshevT(1:N,t1bar(1))*I(2:end,3));
    r0/W(2) - I_C +  I(1,3) + 2*ChebyshevT(1:N,t1bar(2))*I(2:end,3)];
check = min(certificate<0);
% Build the enclosures
r0t = r0/W(3);
w = intval(para(6)); 
p = intval(p);
t0_left = (2*(p(1)-r0t) + w*(t0bar(1)-1))/2;
t0_right = (2*(p(1)+r0t) + w*(t0bar(2)-1))/2;
t0 = infsup(inf(t0_left),sup(t0_right));
t1_left = (p(1)+p(2) - 2*(1+abs(t1bar(1)))*r0t + t1bar(1)*(p(2)-w-p(1)) - w)/2;
t1_right = (p(1)+p(2) + 2*(1+abs(t1bar(2)))*r0t + t1bar(2)*(p(2)-w-p(1)) - w)/2;
t1 = infsup(inf(t1_left),sup(t1_right));
end
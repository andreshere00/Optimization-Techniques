function [a,b,k,sol] = bisection_method(f,a,b,e)
% bisection search method

k = 0;                  % number of iterations
max_iter = 50;          % maximum number of iterations
sol = [];

syms x;
df = diff(f,x);         % computing derivative

while ((abs(b-a)>e) && (k<max_iter))
    k= k + 1;
    c=(a + b)/2;
    dfc = subs(df,x,c);
    dfa = subs(df,x,a);
    if sign(dfc) == sign(dfa)
        a=c;
    else
        b=c;
    end
    sol = [sol, c];
    if k == max_iter
        error('The algorithn could not converge.')
    end
end
end
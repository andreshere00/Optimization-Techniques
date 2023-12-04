function [a,b,k,sol] = fibonacci_method(f,a,b,e)
% Fibonacci search method

k = 0;          % number of iterations 
n = (b-a)/e;    % criteria to obtain the last fibonacci term   
F = [1,2];      % initial values of fibonacci series

while F(end) < n
    F = [F,F(end)+F(end-1)];
end

sol = [];
for i = 1:length(F)-1 
    k = k+1;
    rho = 1 - F(end-i)/F(end-i+1);
    x_1 = a + rho * (b - a);
    x_2 = a + (1 - rho) * (b - a);
    if (f(x_1) > f(x_2))
        a = x_1;
        sol = [sol,a];
    else
        b = x_2;
        sol = [sol,b];
    end
end
end
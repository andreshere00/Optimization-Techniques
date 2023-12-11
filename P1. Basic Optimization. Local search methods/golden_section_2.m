function [a,b,k,sol] = golden_section_2(f,a,b,e)
% Golden section search method

k = 0;                  % number of iterations
max_iter = 50;          % maximum number of iterations
gr = (sqrt(5)-1)/2;     % golden ratio

sol = [];
iter = ceil((log(a-b)-log(e))/log((1+sqrt(5))/2));

while (k < iter)
    k = k + 1;
    x_1 = a + gr*(b-a);
    x_2 = a + (1-gr)*(b-a);
    if (f(x_1) > f(x_2))
        a = x_1;
        sol = [sol,a];
    else
        b = x_2;
        sol = [sol,b];
    end
end
end
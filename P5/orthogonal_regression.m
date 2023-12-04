function output = orthogonal_regression(x,y,E)
% Orthogonal regression optimization problem
    aux = 0;
    for i=1:length(x)
        aux = aux + ((E(1) + E(2).*x(i) + E(3).*x(i).^2 - y(i)).^2 ...
            + (x(i) - E(i+3)).^2 + (x(i).^2 - E(i+3).^2).^2);
    end
    output = aux;
end
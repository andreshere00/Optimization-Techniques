function result = primaldual(A, b, c, eps, rho, sigma, niter)
%% Primal-dual interior-point algorithm
x=A'*inv(A*A')*b;
u=inv(A*A')*A*c;
s=c-A'*u;

minx=min(x);
mins=min(s);
delx=max(-1.5*minx,0);
dels=max(-1.5*mins,0);

x=x+delx;
s=s+dels;
delxx=(1/2)*(x'*s)/sum(s);
delss=(1/2)*(x'*s)/sum(x);

x=x+delxx;
s=s+delss;

step=1;
aux=1;
while(step>eps || aux<=niter)

    n=length(x);
    mu=(sigma*s'*x)/n;
    e=ones(1, n)';
    X=diag(x);
    S=diag(s);
    I=eye(n);

    rows_S = size(S, 1); cols_S = size(S, 2);
    rows_X = size(X, 1); cols_X = size(X, 2);
    rows_A = size(A, 1); cols_A = size(A, 2); cols_A_transpose = size(A', 2);

    result = [S, zeros(rows_S, cols_X), X;
        A, zeros(rows_A, cols_X), zeros(rows_A, cols_A);
        zeros(cols_A_transpose, cols_S), A', I];

    temp = inv(result)*[mu*e-X*S*e;b-A*x;c-A'*u-s];

    dx = temp(1:size(S, 2));
    du = temp(size(S, 2) + 1:size(S, 2) + size(X, 2));
    ds = temp(size(S, 2) + size(X, 2) + 1:end);

    if ds >= 0
        disp('It does not exist a feasible solution');
        return;
    end

    if dx >= 0
        disp('The solution is not bounded');
        return;
    end

    dxneg = dx(dx < 0);
    for i = 1:numel(dxneg)
        alpha_p(i) = rho * min(-x ./ dxneg);
    end
    alpha_p = min(alpha_p);

    dsneg = ds(ds < 0);
    for i = 1:numel(dsneg)
        alpha_d(i) = rho .* min(-s ./ dsneg);
    end
    alpha_d = min(alpha_d);

    s = s + alpha_d.*ds;
    u = u + alpha_d.*du;
    x = x + alpha_p.*dx;

    step=transpose(c)*x-transpose(b)*u;
    aux = aux+1;
    result = x;
end
end
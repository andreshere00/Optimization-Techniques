function [sol,flag, prompt] = primal_dual_interior_point(A, b, c, eps, rho, sigma, niter)
%% primal_dual_interior_point: Solves linear programming problems using the primal-dual interior-point algorithm.
%
%   [sol, flag, prompt] = PRIMAL_DUAL_INTERIOR_POINT(A, b, c, eps, rho, sigma, niter) solves the linear programming problem
%
%      minimize  c' * x
%      subject to  A * x = b
%                  x >= 0
%
%   using the primal-dual interior-point algorithm. It returns the optimal solution 'sol', a flag 'flag', and a status message 'prompt'.
%
%   Input:
%      A      - Coefficient matrix for equality constraints (m x n)
%      b      - Right-hand side of equality constraints (m x 1)
%      c      - Coefficient vector for the objective function (n x 1)
%      eps    - Tolerance for stopping criterion
%      rho    - Step size parameter
%      sigma  - Barrier parameter
%      niter  - Maximum number of iterations
%
%   Output:
%      sol    - Optimal solution to the linear program (n x 1)
%      flag   - Flag indicating the result:
%               -2: Infeasible solution
%               -1: Unbounded solution
%                1: Optimal solution found
%      prompt - Status message describing the result

% Initialization
x0=A.'*inv(A*A.')*b;
u0=inv(A*A.')*A*c;
s0=c-A.'*u0;

% Minimum
minx=min(x0);
mins=min(s0);
delx=max(-(3/2)*minx,0);
dels=max(-(3/2)*mins,0);

x0=x0+delx;
s0=s0+dels;

delxx=(1/2)*(x0.'*s0)/(ones(size(x0)).'*s0);
delss=(1/2)*(x0.'*s0)/(ones(size(s0)).'*x0);

x0=x0+delxx;
s0=s0+delss;

% Initial solution, first step
x=x0;
u=u0;
s=s0;

% Iterations
step=Inf;
iter=0;
while (step > eps) | (iter < niter)

    % Initial conditions of each step.
    iter = iter+1;
    n = length(x);
    m = length(u);
    X = diag(x);
    S = diag(s);
    mu = sigma*(s.'*x)/n;

    % dw = A*Z; with A = inv(M).
    bb1=mu*ones(n,1)-X*S*ones(n,1);
    bb2=b-A*x;
    bb3=c-A.'*u-s;

    M=[S, zeros(n,m), X; A, zeros(m,m), zeros(m,n); zeros(n,n), A.', eye(n)];

    Z=[bb1;bb2;bb3];

    dw = inv(M)*Z;

    % Subvectors which will evaluate each step of the algorithm.
    % dw = [dx, du, ds]
    dx=zeros(1,n);
    for i=1:n
        dx(i)=dw(i);
    end
    dx=dx.';

    du=zeros(1,n+m);
    for i=(n+1):(n+m)
        du(i)=dw(i);
    end
    du(1:n)=[]; du=du.';

    ds=zeros(1,n+m+n);
    for i=(n+m+1):(n+m+n)
        ds(i)=dw(i);
    end
    ds(1:(n+m))=[]; ds=ds.';

    % Evaluation
    if all(ds>=0)
        prompt = 'Feasible solution does not exists. Initial point returned.';
        flag = -2;
        sol = x0;
        break
    elseif all(dx >= 0)
        prompt = 'Not bounded solution. -Inf returned.';
        flag = -1;
        sol = -Inf;
        break
    else

        % Computing alpha_p and alpha_d 
        alpha_p_vector=-x(dx < 0)./dx(dx < 0);
        alpha_d_vector=-s(ds < 0)./ds(ds < 0);

        alpha_p=rho * min(alpha_p_vector); % Primal step magnitude.
        alpha_d=rho * min(alpha_d_vector); % Dual step magnitude.

        % New slack variable for the dual problem
        s = s + alpha_d.*ds;
        % New solution for the dual problem
        u = u + alpha_d.*du;
        % New solution for the primal problem
        x = x + alpha_p.*dx;
        % step
        step = (c.'*x-b.'*u);
        prompt = sprintf('Iteration %d; Duality step: %.4f.', iter, step);
        flag = 1;
        sol = x.';
    end
end
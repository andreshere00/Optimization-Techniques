% Define the linear programming problem:
% Minimize: c'*x
% Subject to: A*x = b, x >= 0

A=[-1,2;2,1];
b=[4;6];
c=[-1;-1];

% initializing variables
x=transpose(A)*inv(A*transpose(A))*b;
u=inv(A*transpose(A))*A*c;
s=c-transpose(A)*u;

n = length(x);
I = eye(n);

minx=min(x);
mins=min(s);
delx=max(-1.5*minx,0);
dels=max(-1.5*mins,0);
x=x+delx
s=s+dels
delxx=(0.5)*(x'*s)/s;
delss=(0.5)*(x'*s)/x;
x=x+delxx
s=s+delss

%% first iteration
sigma = 0.5;
rho = -0.5;
eps = 0.5;

mu=(sigma*s'*x)/n;
e=ones(n,1);
X=diag(x);
S=diag(s);
I=eye(n);

bb1 = mu.*e - X.*S.*e; 
bb2 = b - A*x;
bb3 = c - transpose(A).*u -s;
bb = [bb1;bb2;bb3];

M1 = [S, zeros(2,4), X];
M2 = [A, zeros(size(bb2,2),size(bb2,1)), zeros(size(bb3,1),size(bb3,2))];
M3 = [zeros(size(bb2,2),size(bb2,1)),transpose(A),I];

M = [M1;M2;M3];

inv(M)
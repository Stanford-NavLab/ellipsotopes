%% description

%% parameters

% dimension
n = 10;

% iterations
N = 10;

% some linear system
A = rand(n,n);

x0 = zeros(n,1);
cov = 0.1*eye(n);

% zonotope reachability
x0_zono = zonotope(x0,cov);
x_zono = x0_zono;

for i = 1:N
    x_zono = A * x_zono;
end
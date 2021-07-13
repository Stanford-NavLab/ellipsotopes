%% description
% Fault detection example
%
% References:
% [1] Scott, J.K. Constrained zonotopes
%
%% setup
% linear system
A = [2.0 1.4; 
     0.5 3.1];
B = eye(2);

C = [0.2 1.0;
     0.8 0.3];
D = eye(2);

x0 = [0;0];
N = 100; % number of iterations

% propagate system for N iterations and generate (noiseless) measurements
x = x0;
Y = zeros(2,N);
for i = 1:N
    x = A*x;
    y = C*x;
    Y(:,i) = y;
end

% initial set of states
X0 = ellipsotope(2,[0;0],0.1*eye(2));

% noise sets
W = ellipsotope(2,[0;0],0.1*eye(2));
V = ellipsotope(2,[0;0],0.1*eye(2));

% set-based estimator from [1] eqn (32)
Xi = X0; X = cell(1,N);
for i = 1:N
    Xi = intersect(A*Xi + B*W, Y(:,i) + (-1)*D*V, C);
    X{i} = Xi;
end

figure(1); hold on
for i = 1:N
    plot(X{i});
end

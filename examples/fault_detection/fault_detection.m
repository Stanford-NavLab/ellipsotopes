%% description
% Fault detection example

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
X = ellipsotope(2,[0;0],0.1*eye(2));
% noise sets
W = ellipsotope(2,[0;0],0.1*eye(2));
V = ellipsotope(2,[0;0],0.1*eye(2));

% propagate set-based system while performing fault check
for i = 1:N
    X = (A*X + B*W) & (Y(:,i) + (-1)*D*V);
end

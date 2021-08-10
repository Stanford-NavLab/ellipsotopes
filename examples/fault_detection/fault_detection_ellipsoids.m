%% description
% Fault detection example using ellipsoids
%
% References:
% [1] Scott, J.K. Constrained zonotopes
%
clear 
%% setup
% Nominal (1) and Faulty (2) model parameters
Ra(1) = 1.2030; Ra(2) = 1.5030;
L(1) = 5.5840; L(2) = 5.5840;
Ke(1) = 8.5740; Ke(2) = 8.5740;
Kt = 1.0005 * Ke;
J1(1) = 1.4166; J1(2) = 1.4166;
fr(1) = 2.4500; fr(2) = 2.4500;

% time discretization
dt = 1e-3;

% Nominal model
A{1} = dt * [-Ra(1)/L(1) -Ke(1)/L(1); 
             Kt(1)/J1(1) -fr(1)/J1(1)];
B{1} = dt * [1/L(1); 0];
Bw{1} = [-0.0085 -0.0006;
         -0.0603 0.0002];

% Faulty model
A{2} = dt * [-Ra(2)/L(2) -Ke(2)/L(2); 
             Kt(2)/J1(2) -fr(2)/J1(2)];
B{2} = dt * [1/L(2); 0];
Bw{2} = [-0.0101 -0.0006;
         -0.0595 0.0002];

C = eye(2);
D = eye(2);

% Nominal steady-state
x_N = [0.2; 70.3]; % A, rad/s
u_N = 6; % V
% control input limits
u_min = 0; u_max = 12;
% LQ gain
K{1} = dlqr(A{1},B{1},eye(2),0.1);
K{2} = dlqr(A{2},B{2},eye(2),0.1);

% noise sets
W = ellipsoid(eye(2),[0;0]);
V = ellipsoid([0.06 0; 0 0.6],[0;0]);

% initial set of states
X0 = ellipsoid([0.06 0; 0 0.6],[0.6;70]);

%% run simulations

N_sims = 10;
fault_steps = zeros(N_sims,1);
avg_detect_steps = 0;

for i = 1:N_sims
    % sample initial state
    x_0 = sample(X0,1);
    % initial measurement
    % sample v_0
    v_0 = sample(V,1);
    y_0 = C * x_0 + D * v_0;
    N = 100; % number of iterations

    x_k = x_0; y_k = y_0;
%     X = zeros(2,N);
%     Y = zeros(2,N);

    % set-based estimator from [1] eqn (32)
    % initialization
    O_k = X0 & (y_0 + (-1)*D*V); 
    O = cell(1,N);
    fault = 0;

    % apply estimator to nominal model and simulate faulty model 
    for k = 1:N
        % simulate faulty model
        % sample w_k and v_k from W and V
        w_k = sample(W,1);
        v_k = sample(V,1);
        % control law
        u_k = u_N - K{2} * (y_k - x_N);
        % apply saturation limits
        u_k = max(min(u_max, u_k), u_min);
        % state update
        x_k = A{2} * x_k + B{2} * u_k + Bw{2} * w_k;
        % measurement
        y_k = C * x_k + D * v_k;
        % record state and measurement
%         X(:,k) = x_k;
%         Y(:,k) = y_k;
        
        tic
        % fault detection step
        F = C * (A{1}*O_k + Bw{1}*W) + D*V;
        if ~F.containsPoint(y_k)
            fault = k;
            fault_steps(i) = fault;
            break
        end

        % set-based estimator update
        O_k = (A{1}*O_k + Bw{1}*W) & (y_k + (-1)*D*V);
        O{k} = O_k;
        disp(toc)
    end

    avg_detect_steps = avg_detect_steps + fault;
end

disp(['Average timesteps for detection: ', num2str(avg_detect_steps/N_sims)])
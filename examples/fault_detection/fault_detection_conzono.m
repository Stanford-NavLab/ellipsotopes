%% description
% Fault detection example using constrained zonotopes
%
% References:
% [1] Scott, J.K. Constrained zonotopes
%
clear; clc; close all
%% for debugging
%rng(3)

plot_flag = false;

%% setup

load('samples.mat')

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
% A{1} = dt * [-Ra(1)/L(1) -Ke(1)/L(1); 
%              Kt(1)/J1(1) -fr(1)/J1(1)];
% B{1} = dt * [1/L(1); 0];
Bw{1} = [-0.1 -0.2;
         -0.2 0.1];
A{1} = dt * eye(2);
B{1} = dt * [1; 0];

% Faulty model
% A{2} = dt * [-Ra(2)/L(2) -Ke(2)/L(2); 
%              Kt(2)/J1(2) -fr(2)/J1(2)];
% B{2} = dt * [1/L(2); 0];
Bw{2} = [-0.2 -0.2;
         -0.1 0.1];
A{2} = dt * 2 * eye(2);
B{2} = dt * [2; 0];

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
% zonotope noise
%W = conZonotope([0;0],eye(2));
%V = conZonotope([0;0],[0.06 0; 0 0.6]);
% ellipsoidal noise
W_e = ellipsoid(eye(2),[0;0]);
V_e = ellipsoid([0.06 0; 0 0.6],[0;0]);
W = zonotope(W_e,14,'o:norm');
V = zonotope(V_e,14,'o:norm');


% initial set of states
X0 = conZonotope([0.6;70],[0.06 0; 0 0.6]);

% desired order reduction values
n_c = 3; % num constraints
o_d = 5; % degrees-of-freedom order

N_sims = 5; % number of simulations to run
N = 100; % number of iterations

%% run simulations

fault_steps = zeros(N_sims,1);
avg_detect_steps = 0;
avg_step_time = zeros(N_sims,1);
missed_detections = 0;

if plot_flag
    % state domain plot
    f1 = figure(1); axis equal; hold on
    % measurement domain plot
    f2 = figure(2); axis equal; hold on
end

for i = 1:N_sims
    disp(['i = ',num2str(i)])
    % sample initial state
    x_0 = randPoint(X0);
    %x_0 = samples.x0(:,i);
    % initial measurement
    % sample v_0
    v_0 = randPoint(V_e);
    %v_0 = samples.v(i,:,1)';
    y_0 = C * x_0 + D * v_0;

    x_k = x_0; y_k = y_0;

    % set-based estimator from [1] eqn (32)
    % initialization
    O_k = (C * X0) & (y_0 + (-1)*D*V); 
    O = cell(1,N);
    fault = 0;

    % apply estimator to nominal model and simulate faulty model 
    for k = 1:N
        disp([' k = ',num2str(k)])
        % simulate faulty model
        % sample w_k and v_k from W and V
        w_k = randPoint(W_e);
        v_k = randPoint(V_e);
        %w_k = samples.w(i,:,k)';
        %v_k = samples.v(i,:,k+1)';
        disp(['  w_k: (',num2str(w_k(1)),', ',num2str(w_k(2)),') v_k: (',num2str(v_k(1)),', ',num2str(v_k(2)),')']);
        % control law
        u_k = u_N - K{2} * (y_k - x_N);
        % apply saturation limits
        u_k = max(min(u_max, u_k), u_min);
        % state update
        x_k = A{2} * x_k + B{2} * u_k + Bw{2} * w_k;
        % measurement
        y_k = C * x_k + D * v_k;
        disp(['  x_k: (',num2str(x_k(1)),', ',num2str(x_k(2)),') y_k: (',num2str(y_k(1)),', ',num2str(y_k(2)),')']);
        
        if plot_flag
            figure(1); scatter(x_k(1),x_k(2));
            figure(2); scatter(y_k(1),y_k(2));
        end
        
        tic
        % fault detection step
        F = C * (A{1}*O_k + Bw{1}*W) + D*V;
        if ~in(F,y_k)
            disp('Fault detected');
            fault = k;
            fault_steps(i) = fault;
            break
        end

        % set-based estimator update
        O_k = (C * (A{1}*O_k + Bw{1}*W)) & (y_k + (-1)*D*V);
        O{k} = O_k;
        
        % order reduction 
        O_k = reduce(O_k,'scott',o_d); % reduce to o_d degrees of freedom order
        O_k = reduceConstraints(O_k,n_c); % reduce to n_c constraints
        disp([' n_c: ',num2str(size(O_k.A,1)),' n_g: ',num2str(size(O_k.A,2))])

        if plot_flag
            figure(1); plot(O_k);
            figure(2); plot(F);
        end
        
        %disp(toc)
        avg_step_time(i) = avg_step_time(i) + toc;
        
        if plot_flag
            clf(f1); clf(f2); % set breakpoint here for plotting
        end
    end
    
    if fault == 0
        disp('Failed to detect fault');
        missed_detections = missed_detections + 1;
    else
        avg_detect_steps = avg_detect_steps + fault;
    end
    avg_step_time(i) = avg_step_time(i) / k;
end

disp(['Average timesteps for detection: ', num2str(avg_detect_steps/(N_sims-missed_detections))])
disp(['Average time per timestep: ', num2str(mean(avg_step_time))])
disp(['Missed detections: ', num2str(missed_detections)])
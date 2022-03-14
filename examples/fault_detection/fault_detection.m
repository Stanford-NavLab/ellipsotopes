%% description
% Fault detection example
%
% References:
% [1] Scott, J.K. Constrained zonotopes
%
% Authors: Adam Dai
% Created: Shrug
% Updated: 14 Mar 2022 (Shreyas did some debuggins)
clear ; clc ;
%% user parameters
% random number generator seed (for replicability)
rng(0)

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
% zonotope noise
%W = ellipsotope(2,[0;0],eye(2),[],[],{1,2});
%V = ellipsotope(2,[0;0],[0.06 0; 0 0.6],[],[],{1,2});
% ellipsoidal noise
W = ellipsotope(2,[0;0],eye(2));
V = ellipsotope(2,[0;0],[0.06 0; 0 0.6]);

% initial set of states
X0 = ellipsotope(2,[0.6;70],[0.06 0; 0 0.6],[],[],{1,2});

% desired order reduction values
n_c = 3; % num constraints
o_d = 5; % degrees-of-freedom order
n_g = 13; % corresponding number of generators (for dimension 2)

% number of simulations to run
N_sims = 10 ;

%% run simulations
% set up to save info
fault_steps = zeros(N_sims,1);
avg_detect_steps = 0;
avg_step_time = zeros(N_sims,1);
missed_detections = 0;

for i = 1:N_sims
    % sample initial state
    x_0 = sample_from_ellipsotope(X0);
    % initial measurement
    % sample v_0
    v_0 = sample_from_ellipsotope(V);
    y_0 = C * x_0 + D * v_0;
    N = 100; % number of iterations

    x_k = x_0; y_k = y_0;
%     X = zeros(2,N);
%     Y = zeros(2,N);

    % set-based estimator from [1] eqn (32)
    % initialization
    O_k = intersect(X0, y_0 + (-1)*D*V, C); 
    O = cell(1,N);
    fault = 0;

    % apply estimator to nominal model and simulate faulty model 
    for k = 1:N
        % simulate faulty model
        % sample w_k and v_k from W and V
        w_k = sample_from_ellipsotope(W);
        v_k = sample_from_ellipsotope(V);
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
%         F_cond = cond(F.constraint_A) ;
%         
%         if F_cond > 1e8
%             dbstop in fault_detection at 116
%             disp('hi')
%         end
        
        %%% SHREYAS DEBUGGING DEGENERATE CONSTRAINTS %%%
        
%         A_test = F.constraint_A ;
%         b_test = F.constraint_b ;
%         x_test = pinv(A_test)*b_test ;
%         if vecnorm(A_test*x_test - b_test) > 1e-7
%             dbstop in fault_detection at 120
%             disp('hi')
%         end
        %%% SHREYAS DEBUGGING DEGENERATE CONSTRAINTS %%%
        
        if ~F.contains(y_k)
            fault = k;
            fault_steps(i) = fault;
            break
        end

        % set-based estimator update
        O_k = intersect(A{1}*O_k + Bw{1}*W, y_k + (-1)*D*V, C);
        O{k} = O_k;
        
        % order reduction
        %O_k = reduce_2_etope_to_minimal_exact_rep(O_k); 
        O_k = reduce_constraint(O_k,n_c); % reduce to n_c constraints
        disp(['n_c: ',num2str(size(O_k.constraint_A,1)),' n_g: ',num2str(size(O_k.constraint_A,2))])
        
        disp(toc)
        avg_step_time(i) = avg_step_time(i) + toc;
    end
    
    if fault == 0
        disp('Failed to detect fault');
        missed_detections = missed_detections + 1;
    else
        avg_detect_steps = avg_detect_steps + fault;
    end
    avg_step_time(i) = avg_step_time(i) / k;
end

disp(['Average timesteps for detection: ', num2str(avg_detect_steps/N_sims)])
disp(['Average time per timestep: ', num2str(mean(avg_step_time))])
disp(['Missed detections: ', num2str(missed_detections)])
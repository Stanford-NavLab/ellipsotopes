%% description
% Linear reachability
% Given a robot running an LQR controler and KF state estimator, and
% nominal trajectory, compute forward reachable set

clc
clear
% close all

%% setup

% system
dx = 1; dt = 0.2; N = 50;
A = eye(2); B = dt*eye(2); C = eye(2);
K = dlqr(A,B,eye(2),eye(2));
Q = diag([0.05 0.05]);
R = diag([0.05 0.05]);

n = size(A,1); % system dimension
m = size(R,1); % measurement dimension

% initial state
x0 = [0;0]; P0 = 0.1*eye(2);

% nominal trajectory
u_nom = ones(n,N-1);
x_nom = zeros(n,N); x_nom(:,1) = x0;
for i = 2:N
    x_nom(:,i) = A*x_nom(:,i-1) + B*u_nom(:,i-1);
end

N_traj = 100;
%% generate sample trajectories
X = nan(2,N,N_traj);
Z = nan(2,N,N_traj);

for i_traj = 1:N_traj
    
    % initialize state and covariance
    x_est = x0; P = P0;
    x = mvnrnd(x_est', P0, 1)';

    % initialize trajectory
    X(:,1,i_traj) = x;

    for i = 2:N

        % apply feedback control
        err = x_est - x_nom(:,i-1);
        u = u_nom(:,i-1) - K*err;

        % dynamics
        w = mvnrnd(zeros(n,1), Q, 1)';
        x = A*x + B*u + w;
        X(:,i,i_traj) = x;

        % noisy/biased measurement
        v = mvnrnd(zeros(m,1), R, 1)';
        z = C*x + v;

        % Kalman filter predict
        x_pred = A*x_est + B*u;
        P_pred = A*P*A' + Q;

        % Kalman filter update
        L = P_pred*C'/(C*P_pred*C' + R);
        x_est = x_pred + L*(z - C*x_pred);
        P = P_pred - L*C*P_pred;

    end
    
end

%% ellipsotope




%% reachability calculation
P = P0;
pP0 = probZonotope([0;0],cov2probGen(P),3);

WpZ = probZonotope([0;0],cov2probGen(Q),3);
VpZ = probZonotope([0;0],cov2probGen(R),3);

pXrs = cell(1,N);
pXrs{1} = pP0;

coeff_a = eye(2); coeff_b = 0;
coeff_c = cell(1,N); coeff_d = cell(1,N);
coeff_c{1} = nan; coeff_d{1} = nan;

coeff_e = eye(2);
coeff_p = cell(1,N); coeff_q = cell(1,N);
coeff_p{1} = nan; coeff_q{1} = nan;

for k = 2:N
    
    %update coeffs a and b
    coeff_a = (A-B*K)*coeff_a;
    coeff_b = (A-B*K)*coeff_b - B*K*coeff_e;
    
    %update coeffs c and d
    for n = 2:k-1
        coeff_c{n} = (A-B*K)*coeff_c{n} + -B*K*coeff_p{n};
        coeff_d{n} = (A-B*K)*coeff_d{n} + -B*K*coeff_q{n};
    end
    
    %add new coeffs c and d
    coeff_c{k} = eye(2);  coeff_d{k} = zeros(2);
    
    %calculate all CpZ and DpZ terms
    all_CpZ = coeff_c{k}*WpZ;
    all_DpZ = coeff_d{k}*VpZ;
    for n = 2:k-1
        all_CpZ = all_CpZ + coeff_c{n}*WpZ;
        all_DpZ = all_DpZ + coeff_d{n}*VpZ;
    end
    
    %compute reachable set
    pXr = ( coeff_a - coeff_b )*pP0 + all_CpZ + all_DpZ;
%     pXr = reduce(pXr,'girard',pZ_order);
    pXrs{k} = pXr + x_nom(:,k);
    
    %online filter steps
    P_pred = A*P*A' + Q;
    L = P_pred*C'/(C*P_pred*C' + R);
    P = P_pred - L*C*P_pred;
    
    %update coeff e
    coeff_e = (eye(2) - L*C)*A*coeff_e;
    
    %update coeffs p and q
    for n = 2:k-1
        coeff_p{n} = (eye(2) - L*C)*A*coeff_p{n};
        coeff_q{n} = (eye(2) - L*C)*A*coeff_q{n};
    end
    
    %add coeffs p and q for new w and v
    coeff_p{k} = -(eye(2) - L*C);  coeff_q{k} = L;
end

%% plot 3-sigma confidence zonotopes
m = 3;
Xrs = cell(1,N);
figure()
hold on; grid on;
axis equal;
for i_traj = 1:N_traj
    plot(X(1,:,i_traj),X(2,:,i_traj),'b');
end
for i = 1:N
    umeanZ = zonotope(pXrs{i}.Z);
    covZ = cov2zonotope(pXrs{i}.cov,m,2);
    Xrs{i} = umeanZ + covZ;
    plot(Xrs{i});
end
axis equal
xlabel('x-coordinate (m)');
ylabel('y-coordinate (m)');
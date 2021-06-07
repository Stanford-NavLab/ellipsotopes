%% description
% Dubins reachability
% Given a robot running an LQR controler and KF state estimator, and
% nominal trajectory, compute forward reachable set

clc
clear
% close all

%% robot matrices

params.dt = 0.2;
params.V = 5;
params.dx = params.V*params.dt;
%params.turn_radius = params.V/(20*pi/180);
params.turn_radius = 12;
params.z_cov_range = 0.1;
params.z_bias_range = 1;
params.z_cov_theta = 0.001;
params.pZ_order = 10;
params.N_rollouts = 100;
params.m = 3;

% motion and measurement uncertainty
robot.Q_lqr = diag([0.1 0.1 10]);
robot.R_lqr = diag([10 0.1]);
robot.Q = 10*diag([0.01, 0.01, 0.001]);

range_sigma = 0.4;
heading_sigma = 0.001;
robot.R = diag([range_sigma*ones(1,4), heading_sigma]);
robot.P0 = diag([0.1 0.1 0.01]);

% robot dimensions
robot.width = 1;
robot.length = 2;


reach.beacon_positions = [[-10;-10], [60;-10], [60;60], [-10;60]];

%% generate nominal trajectory

trajectory.waypoints = [ [0; 0; pi/2], [40; 10; pi/4], [50; 50; pi/8] ];
% trajectory.waypoints = [ [0; 0; -pi/2], [10; -30; 0] ];
trajectory.dubins_offset = 0;
[trajectory.x_nom, trajectory.u_nom, trajectory.K_nom, ~, trajectory.dubins_offset] = dubins_curve(trajectory.waypoints(:,1)', trajectory.waypoints(:,2)', params.turn_radius, params.dx, params.dt, trajectory.dubins_offset, 1, robot);
for i = 3:size(trajectory.waypoints,2)
    [trajectory.x_nom_t, trajectory.u_nom_t, trajectory.K_nom_t, ~, trajectory.dubins_offset] = dubins_curve(trajectory.x_nom(:,end), trajectory.waypoints(:,i), params.turn_radius, params.dx, params.dt, trajectory.dubins_offset, 1, robot);
    trajectory.x_nom = cat(2, trajectory.x_nom, trajectory.x_nom_t);
    trajectory.u_nom = cat(2, trajectory.u_nom, trajectory.u_nom_t);
    trajectory.K_nom = cat(3, trajectory.K_nom, trajectory.K_nom_t);
end
trajectory.N_timesteps = size(trajectory.x_nom,2);

for k = 2:trajectory.N_timesteps
    trajectory.x_nom(:,k) = trajectory.x_nom(:,k-1) + [trajectory.u_nom(1,k-1)*cos(trajectory.x_nom(3,k-1)); trajectory.u_nom(1,k-1)*sin(trajectory.x_nom(3,k-1)); trajectory.u_nom(2,k-1)]*params.dt;
end

%% generate trajectory rollouts

rollouts.X = nan(3,trajectory.N_timesteps,params.N_rollouts);

for i = 1:params.N_rollouts
    
    rollouts.x_est = trajectory.x_nom(:,1); rollouts.P = robot.P0;
    rollouts.x = mvnrnd(rollouts.x_est', robot.P0, 1)';
    
    rollouts.X(:,1,i) = rollouts.x;
        
    for k = 2:trajectory.N_timesteps
        
        rollouts.err = rollouts.x_est - trajectory.x_nom(:,k-1);
        rollouts.u = trajectory.u_nom(:,k-1) - trajectory.K_nom(:,:,k-1)*rollouts.err;
        
        rollouts.w = mvnrnd([0 0 0], robot.Q, 1)';
        rollouts.x = (rollouts.x + [rollouts.u(1)*cos(rollouts.x(3)); rollouts.u(1)*sin(rollouts.x(3)); rollouts.u(2)]*params.dt) + rollouts.w;
        rollouts.X(:,k,i) = rollouts.x;
        
        %obtain system matrices
        [robot.A, robot.B, robot.C, robot.K] = obtain_system_matrices( k, trajectory, reach, params );
        v = mvnrnd(zeros(5,1),robot.R);
        for n = 1:4
            rollouts.z(n,1) = sqrt((rollouts.x(1)-reach.beacon_positions(1,n))^2 + (rollouts.x(2)-reach.beacon_positions(2,n))^2 ) + v(n);
        end
        rollouts.z(5,1) = rollouts.x(3) + v(5);
%         robot.A = eye(3);
%         robot.A(1,3) = -trajectory.u_nom(1,k-1)*params.dt*sin(trajectory.x_nom(3,k-1));
%         robot.A(2,3) = trajectory.u_nom(1,k-1)*params.dt*cos(trajectory.x_nom(3,k-1));
%         
%         robot.C = nan(5,3); rollouts.z = nan(5,1);
%         for n = 1:4
%             rollouts.d_nom = sqrt((trajectory.x_nom(1,k)-reach.beacon_positions(1,n))^2 + (trajectory.x_nom(2,k)-reach.beacon_positions(2,n))^2 );
%             robot.C(n,:) = [(trajectory.x_nom(1,k)-reach.beacon_positions(1,n))/rollouts.d_nom (trajectory.x_nom(2,k)-reach.beacon_positions(2,n))/rollouts.d_nom 0];
%             rollouts.z(n,1) = sqrt((rollouts.x(1)-reach.beacon_positions(1,n))^2 + (rollouts.x(2)-reach.beacon_positions(2,n))^2 ) + normrnd(0,sqrt(params.z_cov_range));
%         end
%         robot.C(5,:) = [0 0 1];
%         rollouts.z(5,1) = rollouts.x(3) + normrnd(0,sqrt(params.z_cov_theta));
        
        %prediction step
        rollouts.x_pred = (rollouts.x_est + [rollouts.u(1)*cos(rollouts.x_est(3)); rollouts.u(1)*sin(rollouts.x_est(3)); rollouts.u(2)]*params.dt);
        rollouts.P_pred = robot.A*rollouts.P*robot.A' + robot.Q;
        
        rollouts.z_pred = nan(5,1);
        for n = 1:4
            rollouts.z_pred(n,1) = sqrt( (rollouts.x_pred(1)-reach.beacon_positions(1,n))^2 + (rollouts.x_pred(2)-reach.beacon_positions(2,n))^2 );
        end
        rollouts.z_pred(5,1) = rollouts.x_pred(3);
        
        rollouts.L = rollouts.P_pred*robot.C'/(robot.C*rollouts.P_pred*robot.C' + robot.R);
        rollouts.x_est = rollouts.x_pred + rollouts.L*(rollouts.z - rollouts.z_pred);
        rollouts.P = rollouts.P_pred - rollouts.L*robot.C*rollouts.P_pred;
        
    end
    
    
end

%% RRBT uncertainty propagation and reachable set formation

Sigma = robot.P0;
Lambda = zeros(3);
RRBT = cell(1,trajectory.N_timesteps);
RRBT{1} = Sigma;

reach.Xrs = cell(1,trajectory.N_timesteps);
P = 0.99; % probability threshold
eps = -2*log(1-P);
c = trajectory.x_nom(1:2,1);
G = sqrtm(eps*RRBT{1}(1:2,1:2));
Q = eps*RRBT{1}(1:2,1:2);
reach.Xrs{1} = ellipsotope(2,c,G);

% circumscribing circle
body_radius = norm([robot.width robot.length]) / 2;
robot.circ = ellipsotope(2,zeros(2,1),body_radius*eye(2));

G_robot = rotation_matrix_2D(trajectory.x_nom(3,1)) * 0.5 * diag([robot.length robot.width]);
robot.body = ellipsotope(2,zeros(2,1),G_robot,[],[],{1,2});

reach.Xrs{1} = reach.Xrs{1} + robot.body;

for k = 2:trajectory.N_timesteps
    
    %obtain system matrices
    [robot.A, robot.B, robot.C, robot.K] = obtain_system_matrices( k, trajectory, reach, params );
    
    % covariance prediction
    Sigma = robot.A*Sigma*robot.A' + robot.Q;
    L = Sigma*robot.C'/(robot.C*Sigma*robot.C' + robot.R);
    
    Lambda = (robot.A - robot.B*robot.K)*Lambda*(robot.A - robot.B*robot.K)' + L*robot.C*Sigma;
    Sigma = Sigma - L*robot.C*Sigma;
    
    % RRBT distribution
    RRBT{k} = Sigma + Lambda;
    
    % confidence ellipsotope
    c = trajectory.x_nom(1:2,k);
    G = sqrtm(eps*RRBT{k}(1:2,1:2));
    Q = eps*RRBT{k}(1:2,1:2);
    reach.Xrs{k} = ellipsotope(2,c,G);
    %reach.Xrs{k} = ellipsoid(Q,c);
    
    % robot body
%     G_robot = rotation_matrix_2D(trajectory.x_nom(3,k)) * 0.5 * diag([robot.length robot.width]);
%     robot.body = ellipsotope(2,zeros(2,1),G_robot,[],[],{1,2});
    
    % compute rotated body
    d_theta = erfinv(P) * RRBT{k}(3,3) * sqrt(2);
    robot.rot_body = rotated_body(robot, trajectory.x_nom(3,k), d_theta);
    
    reach.Xrs{k} = reach.Xrs{k} + robot.rot_body;
end

%% obstacle
E_obs = ellipsotope(2,[0;0],5*[1 0.5 -0.5; 0 0.866 0.866],[],[],{1,2,3});
E_obs = E_obs + ellipsotope(2,[0;0],2*diag([1;2]));
E_obs = rotation_matrix_2D(0.6) * E_obs;
E_obs = E_obs + [25;30];

%% plot reachable sets and rollouts

figure(1); hold on; grid on;
for k = 1:4:trajectory.N_timesteps
    reach_h = plot(reach.Xrs{k},'EdgeAlpha',1.0,'FaceColor','b','EdgeColor','b','LineWidth',1.5);
end
for i = 1:params.N_rollouts
    roll_h = plot(rollouts.X(1,:,i),rollouts.X(2,:,i));
    roll_h.Color=[0,0,0,0.5];
end
obs_h = plot(E_obs,'FaceColor','r','EdgeColor','r','FaceAlpha',0.5,'EdgeAlpha',1.0);
beac_h = scatter(reach.beacon_positions(1,:),reach.beacon_positions(2,:),100,'black','^','filled');
axis equal
ax = gca; ax.YAxis.FontSize = 8; ax.XAxis.FontSize = 8;
xlabel('x [m]','Interpreter','latex','FontSize',12);
ylabel('y [m]','Interpreter','latex','FontSize',12);
xlim([min(reach.beacon_positions(1,:))-10 max(reach.beacon_positions(1,:))+10]);
ylim([min(reach.beacon_positions(2,:))-10 max(reach.beacon_positions(2,:))+10]);
legend([reach_h roll_h obs_h beac_h],'Reachable sets','Trajectory rollouts','Obstacle','Ranging beacons');
set(gca,'fontsize',15)
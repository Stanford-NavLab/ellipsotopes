%% description
% Dubins reachability
% Given a robot running an LQR controler and KF state estimator, and
% nominal trajectory, compute forward reachable set

clc
clear
% close all

%% robot matrices

params.dt = 0.1;
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
robot.Q = diag([0.01, 0.01, 0.001]);

robot.P0 = diag([0.1 0.1 0.01]);

% robot dimensions
robot.width = 2;
robot.length = 4;

% measurements
reach.beacon_positions = [[-10;-10], [60;-10], [60;60], [-10;60]];
meas_plane = 30;
range_sigma_1 = 0.4;
range_sigma_2 = 10.0;
heading_sigma = 0.001;
robot.R1 = diag([range_sigma_1*ones(1,4), heading_sigma]);
robot.R2 = diag([range_sigma_2*ones(1,4), heading_sigma]);

%% obstacle
n_obs = 3;

obs.etope{1} = ellipsotope(2,[26;28],8*eye(2),[],[],{1,2});
obs.etope{1} = obs.etope{1} + ellipsotope(2,[0;0],eye(2));
obs.etope{2} = ellipsotope(2,[14;-3],5*eye(2),[],[],{1,2});
obs.etope{2} = obs.etope{2} + ellipsotope(2,[0;0],eye(2));
obs.etope{3} = ellipsotope(2,[60;35],5*eye(2),[],[],{1,2});
obs.etope{3} = obs.etope{3} + ellipsotope(2,[0;0],eye(2));

obs.zono{1} = zonotope([26;28],8*eye(2));
obs.zono{2} = zonotope([14;-3],5*eye(2));
obs.zono{3} = zonotope([60;35],5*eye(2));

obs.ellip{1} = ellipsoid(8*eye(2),[26;28]);
obs.ellip{2} = ellipsoid(5*eye(2),[14;-3]);
obs.ellip{3} = ellipsoid(5*eye(2),[60;35]);

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
        
        % compute measurements
        if rollouts.x(1) < meas_plane 
            R = robot.R1;
        else
            R = robot.R2;
        end
        v = mvnrnd(zeros(5,1),R);
        for n = 1:4
            rollouts.z(n,1) = sqrt((rollouts.x(1)-reach.beacon_positions(1,n))^2 + (rollouts.x(2)-reach.beacon_positions(2,n))^2 ) + v(n);
        end
        rollouts.z(5,1) = rollouts.x(3) + v(5);
        
        %prediction step
        rollouts.x_pred = (rollouts.x_est + [rollouts.u(1)*cos(rollouts.x_est(3)); rollouts.u(1)*sin(rollouts.x_est(3)); rollouts.u(2)]*params.dt);
        rollouts.P_pred = robot.A*rollouts.P*robot.A' + robot.Q;
        
        rollouts.z_pred = nan(5,1);
        for n = 1:4
            rollouts.z_pred(n,1) = sqrt( (rollouts.x_pred(1)-reach.beacon_positions(1,n))^2 + (rollouts.x_pred(2)-reach.beacon_positions(2,n))^2 );
        end
        rollouts.z_pred(5,1) = rollouts.x_pred(3);
        
        rollouts.L = rollouts.P_pred*robot.C'/(robot.C*rollouts.P_pred*robot.C' + R);
        rollouts.x_est = rollouts.x_pred + rollouts.L*(rollouts.z - rollouts.z_pred);
        rollouts.P = rollouts.P_pred - rollouts.L*robot.C*rollouts.P_pred;
    end
    
end

%% RRBT uncertainty propagation and reachable set formation

Sigma = robot.P0;
Lambda = zeros(3);
RRBT = cell(1,trajectory.N_timesteps);
RRBT{1} = Sigma;

reach.zono = cell(1,trajectory.N_timesteps);
reach.ellip = cell(1,trajectory.N_timesteps);
reach.etope = cell(1,trajectory.N_timesteps);

P = 0.9973; % probability threshold
eps = -2*log(1-P);
% initial zonotope
c = trajectory.x_nom(1:2,1);
reach.zono{1} = cov2zonotope(RRBT{1}(1:2,1:2),3,2) + c;
% initial ellipsoid
Q = eps*RRBT{1}(1:2,1:2);
reach.ellip{1} = ellipsoid(Q,c);
% initial ellipsotope
G = sqrtm(eps*RRBT{1}(1:2,1:2));
reach.etope{1} = ellipsotope(2,c,G);

% circumscribing circle
body_radius = norm([robot.width robot.length]) / 2;
robot.circ = ellipsotope(2,zeros(2,1),body_radius*eye(2));
d_theta = erfinv(P) * RRBT{1}(3,3) * sqrt(2);
robot.rot_body = rotated_body(robot, trajectory.x_nom(3,1), d_theta);
robot.rot_body_zono = rotated_body_zono(robot, trajectory.x_nom(3,1), d_theta);
robot.rot_body_ellip = rotated_body_ellip(robot);

reach.etope{1} = reach.etope{1} + robot.rot_body;
reach.zono{1} = reach.zono{1} + robot.rot_body_zono;
reach.ellip{1} = reach.ellip{1} + robot.rot_body_ellip;

tic
for k = 2:trajectory.N_timesteps
    
    %obtain system matrices
    [robot.A, robot.B, robot.C, robot.K] = obtain_system_matrices( k, trajectory, reach, params );
    
    if trajectory.x_nom(1,k) < meas_plane 
        R = robot.R1;
    else
        R = robot.R2;
    end
    
    % covariance prediction
    Sigma = robot.A*Sigma*robot.A' + robot.Q;
    L = Sigma*robot.C'/(robot.C*Sigma*robot.C' + R);
    
    Lambda = (robot.A - robot.B*robot.K)*Lambda*(robot.A - robot.B*robot.K)' + L*robot.C*Sigma;
    Sigma = Sigma - L*robot.C*Sigma;
    
    % RRBT distribution
    RRBT{k} = Sigma + Lambda;
    
    % zonotope position uncertainty bound
    c = trajectory.x_nom(1:2,k);
    reach.zono{k} = cov2zonotope(RRBT{k}(1:2,1:2),3,2) + c;
    
    % ellipsoid position uncertainty bound
    Q = eps*RRBT{k}(1:2,1:2);
    reach.ellip{k} = ellipsoid(Q,c);
    
    % ellipsotope position uncertainty bound
    G = sqrtm(eps*RRBT{k}(1:2,1:2));
    reach.etope{k} = ellipsotope(2,c,G);
    
    % compute rotated body overapproximation
    d_theta = erfinv(P) * RRBT{k}(3,3) * sqrt(2);
    
    % - zonotope
    robot.rot_body_zono = rotated_body_zono(robot, trajectory.x_nom(3,k), d_theta);
    reach.zono{k} = reach.zono{k} + robot.rot_body_zono;
    
    % - ellipsoid
    robot.rot_body_ellip = rotated_body_ellip(robot);
    reach.ellip{k} = reach.ellip{k} + robot.rot_body_ellip;
    
    % - ellipsotope
    robot.rot_body = rotated_body(robot, trajectory.x_nom(3,k), d_theta);
    reach.etope{k} = reach.etope{k} + robot.rot_body;
end
disp(['time to compute reachable set: ',num2str(toc)]);

%% collision check

% % zonotope
% disp('Zonotope collision check');
% tic
% in_collision = false;
% for k = 1:trajectory.N_timesteps
%     for j = 1:n_obs
%         if ~isempty(reach.zono{k} & obs.zono{j})
%             in_collision = true;
%             break
%         end
%     end
% end 
% if in_collision
%     disp('  Reachable set collides with obstacle');
% else
%     disp('  Reachable set does not collide with obstacle');
% end
% disp(['  time to collision check: ',num2str(toc)]);
% 
% % ellipsoid
% disp('Ellipsoid collision check');
% tic
% in_collision = false;
% for k = 1:trajectory.N_timesteps
%     for j = 1:n_obs
%         if ~isempty(reach.ellip{k} & obs.ellip{j})
%             in_collision = true;
%             break
%         end
%     end
% end 
% if in_collision
%     disp('  Reachable set collides with obstacle');
% else
%     disp('  Reachable set does not collide with obstacle');
% end
% disp(['  time to collision check: ',num2str(toc)]);
% 
% % ellipsotope
% disp('Ellipsotope collision check');
% tic
% in_collision = false;
% for k = 1:trajectory.N_timesteps
%     for j = 1:n_obs
%         if ~isempty(reach.etope{k} & obs.etope{j})
%             in_collision = true;
%             break
%         end
%     end
% end 
% if in_collision
%     disp('  Reachable set collides with obstacle');
% else
%     disp('  Reachable set does not collide with obstacle');
% end
% disp(['  time to collision check: ',num2str(toc)]);
% 
% %% compute reach set volumes
% reach.zono_vol = 0;
% reach.ellip_vol = 0;
% reach.etope_vol = 0;
% 
% for k = 1:trajectory.N_timesteps
%    reach.zono_vol = reach.zono_vol + volume(reach.zono{k});
%    reach.ellip_vol = reach.ellip_vol + volume(reach.ellip{k});
%    reach.etope_vol = reach.etope_vol + area(reach.etope{k});
% end

%% plot reachable sets and rollouts

c1 = [156, 52, 235] / 255;
c2 = [52, 64, 235] / 255;
c3 = [3, 186, 252] / 255;

figure(1); hold on; grid on; axis equal
% reachable sets
for k = 1:8:trajectory.N_timesteps
    zono_h = plot(reach.zono{k},[1,2],'Filled',true,'FaceAlpha',0.1,'FaceColor',c1,'EdgeColor',c1,'LineWidth',1.5);
    ellip_h = plot(reach.ellip{k},[1,2],'Filled',true,'FaceAlpha',0.1,'FaceColor',c2,'EdgeColor',c2,'LineWidth',1.5);
    etope_h = plot(reach.etope{k},'EdgeAlpha',1.0,'FaceColor',c3,'EdgeColor',c3,'LineWidth',1.5);
end
% obstacles
for i = 1:n_obs
    obs_h = plot(obs.etope{i},'FaceColor','r','EdgeColor','r','FaceAlpha',0.5,'EdgeAlpha',1.0);
end
% beacons
beac_h = scatter(reach.beacon_positions(1,:),reach.beacon_positions(2,:),100,'black','^','filled');
% measurement region
line([meas_plane,meas_plane],[-30,80],'LineStyle','--');
patch([30,80,80,30],[80,80,-30,-30],'r','FaceAlpha',0.1,'EdgeAlpha',0.0);
axis equal
ax = gca; ax.YAxis.FontSize = 8; ax.XAxis.FontSize = 8;
xlabel('x [m]','Interpreter','latex','FontSize',12);
ylabel('y [m]','Interpreter','latex','FontSize',12);
xlim([min(reach.beacon_positions(1,:))-10 max(reach.beacon_positions(1,:))+10]);
ylim([min(reach.beacon_positions(2,:))-10 max(reach.beacon_positions(2,:))+10]);
legend([zono_h ellip_h etope_h],'Zonotope','Ellipsoid','Ellipsotope');
set(gca,'fontsize',15)
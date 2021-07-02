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
params.N_rollouts = 1000;
params.m = 3;

% motion and measurement uncertainty
robot.Q_lqr = diag([0.1 0.1 10]);
robot.R_lqr = diag([10 0.1]);
robot.Q = 10*diag([0.01, 0.01, 0.001]);

robot.P0 = diag([0.1 0.1 0.01]);

% robot dimensions
robot.width = 1;
robot.length = 2;

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
obs = {};
obs{1} = ellipsotope(2,[26;28],8*eye(2),[],[],{1,2});
obs{1} = obs{1} + ellipsotope(2,[0;0],eye(2));
obs{2} = ellipsotope(2,[14;-3],5*eye(2),[],[],{1,2});
obs{2} = obs{2} + ellipsotope(2,[0;0],eye(2));
obs{3} = ellipsotope(2,[60;35],5*eye(2),[],[],{1,2});
obs{3} = obs{3} + ellipsotope(2,[0;0],eye(2));

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

%% initialize reachability and estimator variables

%motion uncertainty
reach.WpZ = probZonotope([0;0;0],cov2probGen(robot.Q),3);
%sensing uncertainties along trajectory
reach.beacon_positions = [ [-20;10], [-20;-80], [30;10], [30;-80] ];
reach.z_bias_w = repmat([params.z_bias_range;params.z_bias_range;params.z_bias_range;params.z_bias_range;0],1,1,trajectory.N_timesteps);
reach.z_cov = repmat(diag([params.z_cov_range params.z_cov_range params.z_cov_range params.z_cov_range params.z_cov_theta]),1,1,trajectory.N_timesteps);

reach.VpZs = cell(1,trajectory.N_timesteps);
robot.Rhats = zeros(5,5,trajectory.N_timesteps);
for k = 1:trajectory.N_timesteps    
    reach.VpZs{1,k} = probZonotope([[0;0;0;0;0], diag(reach.z_bias_w(:,:,k))],cov2probGen(reach.z_cov(:,:,k)),3);
    robot.Rhats(:,:,k) = diag( (( reach.z_bias_w(:,:,k) + params.m*sqrt( diag(reach.z_cov(:,:,k)) ))/params.m).^2 );
%     robot.Rhats(:,:,k) = diag([((reach.z_bias_w(1,:,k) + params.m*sqrt(reach.z_cov(1,1,k)))/params.m)^2 ((reach.z_bias_w(2,:,k) + params.m*sqrt(reach.z_cov(2,2,k)))/params.m)^2 ((reach.z_bias_w(3,:,k) + params.m*sqrt(reach.z_cov(3,3,k)))/params.m)^2 ((reach.z_bias_w(4,:,k) + params.m*sqrt(reach.z_cov(4,4,k)))/params.m)^2]);
end

reach.pXrs = cell(1,trajectory.N_timesteps);
robot.P = robot.P0;

P = 0.9973; % probability threshold
eps = -2*log(1-P);
c = trajectory.x_nom(1:2,1);
G = sqrtm(eps*RRBT{1}(1:2,1:2));
reach.Xcf{1} = ellipsotope(2,c,G);

% circumscribing circle
body_radius = norm([robot.width robot.length]) / 2;
robot.circ = ellipsotope(2,zeros(2,1),body_radius*eye(2));

G_robot = rotation_matrix_2D(trajectory.x_nom(3,1)) * 0.5 * diag([robot.length robot.width]);
robot.body = ellipsotope(2,zeros(2,1),G_robot,[],[],{1,2});

d_theta = erfinv(P) * RRBT{1}(3,3) * sqrt(2);
robot.rot_body = rotated_body(robot, trajectory.x_nom(3,1), d_theta);

% combine robot body with confidence ellipse
reach.Xrs{1} = reach.Xrs{1} + robot.rot_body;

reach.pXr = probZonotope(trajectory.x_nom(:,1), cov2probGen(robot.P0), 3);
reach.pXrs{1} = reach.pXr;
reach.pXest = probZonotope(trajectory.x_nom(:,1),cov2probGen(zeros(3)),3);
reach.pU = trajectory.u_nom(:,1) + -trajectory.K_nom(:,:,1)*(reach.pXest + -trajectory.x_nom(:,1));

reach.phi_1 = eye(3); reach.phi_2 = 0;
reach.phi_3 = cell(1,0); reach.phi_4 = cell(1,0);
reach.phi_5 = cell(1,0); reach.phi_6 = cell(1,0);
reach.phi_7 = cell(1,0); reach.phi_8 = cell(1,0);

reach.phit_2 = eye(3);
reach.phit_3 = cell(1,0); reach.phit_4 = cell(1,0);
reach.phit_5 = cell(1,0); reach.phit_6 = cell(1,0);
reach.phit_7 = cell(1,0); reach.phit_8 = cell(1,0);

reach.V_sets = cell(1,0);
reach.L_f1 = cell(1,0); reach.L_f2 = cell(1,0);
reach.L_h1 = cell(1,0); reach.L_h2 = cell(1,0);

%% propagate reachability

tic;
for k = 2:trajectory.N_timesteps
    
    %obtain system matrices
    [robot.A, robot.B, robot.C, robot.K] = obtain_system_matrices( k, trajectory, reach, params );
    
    %update coeffs phi_1 and phi_2
    reach.phi_1 = (robot.A-robot.B*robot.K)*reach.phi_1;
    reach.phi_2 = (robot.A-robot.B*robot.K)*reach.phi_2 - robot.B*robot.K*reach.phit_2;
    
    %phi_3: update, remove and add new
    [reach.phi_3, reach.remove_3_idx] = update_pXr_coeffs( reach.phi_3, reach.phit_3, robot.A, robot.B, robot.K, params.norm_thresh_W );
    reach.phi_3(reach.remove_3_idx) = [];
    reach.phi_3{ size(reach.phi_3,2)+1 } = eye(3);
    
    %phi_4: update, remove and add new
    [reach.phi_4, reach.remove_4_idx] = update_pXr_coeffs( reach.phi_4, reach.phit_4, robot.A, robot.B, robot.K, params.norm_thresh_V );
    reach.phi_4(reach.remove_4_idx) = []; 
    reach.V_sets(reach.remove_4_idx) = [];
    reach.phi_4{ size(reach.phi_4,2)+1 } = zeros(3,5); 
    reach.V_sets{ size(reach.V_sets,2)+1 } = reach.VpZs{1,k};

    %phi_5: update, remove and add new
    [reach.phi_5, reach.remove_5_idx] = update_pXr_coeffs( reach.phi_5, reach.phit_5, robot.A, robot.B, robot.K, params.norm_thresh_Lf );
    reach.phi_5(reach.remove_5_idx) = [];
    reach.L_f1(reach.remove_5_idx) = [];
    reach.phi_5{ size(reach.phi_5,2)+1 } = eye(3);
    reach.L_f1{ size(reach.L_f1,2)+1 } = f_lagrange_remainder(reach.pXr, reach.pU, trajectory.x_nom(:,k-1), trajectory.u_nom(:,k-1), params.m, params.dt);
    
    %phi_6: update, remove and add new
    [reach.phi_6, reach.remove_6_idx] = update_pXr_coeffs( reach.phi_6, reach.phit_6, robot.A, robot.B, robot.K, params.norm_thresh_Lf );
    reach.phi_6(reach.remove_6_idx) = [];
    reach.L_f2(reach.remove_6_idx) = [];
    reach.phi_6{ size(reach.phi_6,2)+1 } = zeros(3);
    reach.L_f2{ size(reach.L_f2,2)+1 } = f_lagrange_remainder(reach.pXest, reach.pU, trajectory.x_nom(:,k-1), trajectory.u_nom(:,k-1), params.m, params.dt);
    
    %phi_7: update, remove and add new
    [reach.phi_7, reach.remove_7_idx] = update_pXr_coeffs( reach.phi_7, reach.phit_7, robot.A, robot.B, robot.K, params.norm_thresh_Lh );
    reach.phi_7(reach.remove_7_idx) = [];
    reach.L_h1(reach.remove_7_idx) = [];
    reach.phi_7{ size(reach.phi_7,2)+1 } = zeros(3,5);

    %phi_8: update, remove and add new
    [reach.phi_8, reach.remove_8_idx] = update_pXr_coeffs( reach.phi_8, reach.phit_8, robot.A, robot.B, robot.K, params.norm_thresh_Lh );
    reach.phi_8(reach.remove_8_idx) = [];
    reach.L_h2(reach.remove_8_idx) = [];
    reach.phi_8{ size(reach.phi_8,2)+1 } = zeros(3,5);
    
    %calculate all CpZ, DpZ, LRf1, LRf2, LRh1, LRh2
    reach.all_CpZ = reach.phi_3{end}*reach.WpZ;
    for n = 1:size(reach.phi_3,2)-1
        reach.all_CpZ = reach.all_CpZ + reach.phi_3{n}*reach.WpZ;
    end
    reach.all_DpZ = reach.phi_4{end}*reach.V_sets{end};
    for n = 1:size(reach.phi_4,2)-1
        reach.all_DpZ = reach.all_DpZ + reach.phi_4{n}*reach.V_sets{n};
    end
    reach.all_LRf1 = reach.phi_5{end}*reach.L_f1{end};
    for n = 1:size(reach.phi_5,2)-1
        reach.all_LRf1 = reach.all_LRf1 + reach.phi_5{n}*reach.L_f1{n};
    end
    reach.all_LRf2 = reach.phi_6{end}*reach.L_f2{end};
    for n = 1:size(reach.phi_6,2)-1
        reach.all_LRf2 = reach.all_LRf2 + reach.phi_6{n}*reach.L_f2{n};
    end
    reach.all_LRh1 = zonotope(zeros(3,1));
    for n = 1:size(reach.phi_7,2)-1
        reach.all_LRh1 = reach.all_LRh1 + reach.phi_7{n}*reach.L_h1{n};
    end
    reach.all_LRh2 = zonotope(zeros(3,1));
    for n = 1:size(reach.phi_8,2)-1
        reach.all_LRh2 = reach.all_LRh2 + reach.phi_8{n}*reach.L_h2{n};
    end
    
    %compute reachable set
    reach.pXr = trajectory.x_nom(:,k) + ( reach.phi_1 - reach.phi_2 )*(reach.pXrs{1} + -trajectory.x_nom(:,1)) + reach.all_CpZ + reach.all_DpZ + reach.all_LRf1 + reach.all_LRf2 + reach.all_LRh1 + reach.all_LRh2;
    reach.pXrs{k} = reduce( reach.pXr, 'girard', params.pZ_order );
        
    %online filter steps
    robot.P_pred = robot.A*robot.P*robot.A' + robot.Q;
    robot.L = robot.P_pred*robot.C'/(robot.C*robot.P_pred*robot.C' + robot.Rhats(:,:,k));
    robot.P = robot.P_pred - robot.L*robot.C*robot.P_pred;
    
    %update phit_2
    reach.phit_2 = (eye(3) - robot.L*robot.C)*robot.A*reach.phit_2;
    %phit_3: remove, update, add
    reach.phit_3(reach.remove_3_idx) = [];
    reach.phit_3 = update_pXEr_coeffs( reach.phit_3, robot.A, robot.L, robot.C );
    reach.phit_3{ size(reach.phit_3,2)+1 } = -(eye(3) - robot.L*robot.C);
    %phit_4: remove, update, add
    reach.phit_4(reach.remove_4_idx) = [];
    reach.phit_4 = update_pXEr_coeffs( reach.phit_4, robot.A, robot.L, robot.C );
    reach.phit_4{ size(reach.phit_4,2)+1 } = robot.L;
    %phit_5: remove, update, add
    reach.phit_5(reach.remove_5_idx) = [];
    reach.phit_5 = update_pXEr_coeffs( reach.phit_5, robot.A, robot.L, robot.C );
    reach.phit_5{ size(reach.phit_5,2)+1 } = -(eye(3) - robot.L*robot.C);
    %phit_6: remove, update, add
    reach.phit_6(reach.remove_6_idx) = [];
    reach.phit_6 = update_pXEr_coeffs( reach.phit_6, robot.A, robot.L, robot.C );
    reach.phit_6{ size(reach.phit_6,2)+1 } = (eye(3) - robot.L*robot.C);
    %phit_7: remove, update, add
    reach.phit_7(reach.remove_7_idx) = [];
    reach.phit_7 = update_pXEr_coeffs( reach.phit_7, robot.A, robot.L, robot.C );
    reach.phit_7{ size(reach.phit_7,2)+1 } = robot.L;
    %phit_8: remove, update, add
    reach.phit_8(reach.remove_8_idx) = [];
    reach.phit_8 = update_pXEr_coeffs( reach.phit_8, robot.A, robot.L, robot.C );
    reach.phit_8{ size(reach.phit_8,2)+1 } = -robot.L;
    
    %compute predicted state set
    reach.pXpred = (robot.A-robot.B*robot.K)*reach.pXest + reach.L_f2{end};
    
    %obtain the Lagrange remainder sets
    reach.L_h1{ size(reach.L_h1,2)+1 } = h_lagrange_remainder(reach.pXr, trajectory.x_nom(:,k), reach.beacon_positions, params.m);
    reach.L_h2{ size(reach.L_h2,2)+1 } = h_lagrange_remainder(reach.pXpred, trajectory.x_nom(:,k), reach.beacon_positions, params.m);

    %calculate all CpZ, DpZ, LRf1, LRf2, LRh1, LRh2 terms
    reach.all_CpZ = (reach.phi_3{end} + reach.phit_3{end})*reach.WpZ;
    for n = 1:size(reach.phi_3,2)-1
        reach.all_CpZ = reach.all_CpZ + (reach.phi_3{n} + reach.phit_3{n})*reach.WpZ;
    end
    reach.all_DpZ = (reach.phi_4{end} + reach.phit_4{end})*reach.V_sets{end};
    for n = 1:size(reach.phi_4,2)-1
        reach.all_DpZ = reach.all_DpZ + (reach.phi_4{n} + reach.phit_4{n})*reach.V_sets{n};
    end
    reach.all_LRf1 = (reach.phi_5{end} + reach.phit_5{end})*reach.L_f1{end};
    for n = 1:size(reach.phi_5,2)-1
        reach.all_LRf1 = reach.all_LRf1 + (reach.phi_5{n} + reach.phit_5{n})*reach.L_f1{n};
    end
    reach.all_LRf2 = (reach.phi_6{end} + reach.phit_6{end})*reach.L_f2{end};
    for n = 1:size(reach.phi_6,2)-1
        reach.all_LRf2 = reach.all_LRf2 + (reach.phi_6{n} + reach.phit_6{n})*reach.L_f2{n};
    end
    reach.all_LRh1 = reach.phi_7{end}*reach.L_h1{end};
    for n = 1:size(reach.phi_7,2)-1
        reach.all_LRh1 = reach.all_LRh1 + (reach.phi_7{n} + reach.phit_7{n})*reach.L_h1{n};
    end
    reach.all_LRh2 = reach.phi_8{end}*reach.L_h2{end};
    for n = 1:size(reach.phi_8,2)-1
        reach.all_LRh2 = reach.all_LRh2 + (reach.phi_8{n} + reach.phit_8{n})*reach.L_h2{n};
    end
        
    %compute state estimate set
    reach.pXest = trajectory.x_nom(:,k) + ( reach.phi_1 - reach.phi_2 - reach.phit_2 )*(reach.pXrs{1} + -trajectory.x_nom(:,1)) + reach.all_CpZ + reach.all_DpZ + reach.all_LRf1 + reach.all_LRf2 + reach.all_LRh1 + reach.all_LRh2;    
    
    %compute control input set
    reach.pU = trajectory.u_nom(:,k) + -trajectory.K_nom(:,:,k)*(reach.pXest + -trajectory.x_nom(:,k));
    
%     fprintf([num2str(k) '/' num2str(N) '\n']);
end
disp(['Reach sets computed: ' num2str(toc)]);

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

%% Reachability

Sigma = robot.P0;
Lambda = zeros(3);
RRBT = cell(1,trajectory.N_timesteps);
RRBT{1} = Sigma;

reach.Xrs = cell(1,trajectory.N_timesteps); % reachable set (including robot body)


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
    
    % confidence ellipsotope
    c = trajectory.x_nom(1:2,k);
    G = sqrtm(eps*RRBT{k}(1:2,1:2));
    Q = eps*RRBT{k}(1:2,1:2);
    reach.Xcf{k} = ellipsotope(2,c,G);

    % compute rotated body
    d_theta = erfinv(P) * RRBT{k}(3,3) * sqrt(2);
    robot.rot_body = rotated_body(robot, trajectory.x_nom(3,k), d_theta);
    
    reach.Xrs{k} = reach.Xcf{k} + robot.rot_body;
end
disp(['time to compute reachable set: ',num2str(toc)]);

%% collision check
tic
in_collision = false;
collision = {};
for k = 1:trajectory.N_timesteps
    for j = 1:n_obs
        if ~isempty(reach.Xrs{k} & obs{j})
            in_collision = true;
            break
        end
    end
end 
if in_collision
    disp('Reachable set collides with obstacle');
else
    disp('Reachable set does not collide with obstacle');
end
disp(['time to collision check: ',num2str(toc)]);

%% rollout verification

safe_count = zeros(1,trajectory.N_timesteps);
for i = 1:trajectory.N_timesteps
    for j = 1:params.N_rollouts
        % create robot body at this time instant for this rollout
%         pos = rollouts.X(1:2,i,j);
%         theta = rollouts.X(3,i,j);
%         G_robot = rotation_matrix_2d(theta) * 0.5 * diag([robot.length robot.width]);
%         rollouts.body = ellipsotope(2,pos,G_robot);
        % check if robot body is contained in reachable set
        % for now, just check if center of mass point is within confidence
        % ellipsoid
        pos = rollouts.X(1:2,i,j);
        if reach.Xcf{i}.contains(pos)
            safe_count(i) = safe_count(i) + 1;
        end
    end
end

figure(1); 
plot(safe_count);


%% plot reachable sets and rollouts

figure(2); hold on; grid on;
% reachable sets
for k = 1:8:trajectory.N_timesteps
    reach_h = plot(reach.Xrs{k},'EdgeAlpha',1.0,'FaceColor','b','EdgeColor','b','LineWidth',1.5);
end
% rollouts
for i = 1:params.N_rollouts
    roll_h = plot(rollouts.X(1,:,i),rollouts.X(2,:,i));
    roll_h.Color=[0,0,0,0.5];
end
% obstacle
for i = 1:n_obs
    obs_h = plot(obs{i},'FaceColor','r','EdgeColor','r','FaceAlpha',0.5,'EdgeAlpha',1.0);
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
legend([reach_h roll_h obs_h beac_h],'Reachable sets','Trajectory rollouts','Obstacle','Ranging beacons');
set(gca,'fontsize',15)
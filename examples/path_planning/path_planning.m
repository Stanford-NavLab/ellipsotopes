%% description
% Robot path planning example
% Robot's volume is represented by a zonotopic ellipsotope, and the
% (probabilistically-bounded) uncertain position of the robot is
% represented by an ellipsoidal ellipsotope.
%
% Robot uses Kalman filter to estimate its position. KF state estimate and
% covariance is used to create uncertain position ellipse.

clear;
%% parameters

% sim
dt = 0.1;
goal_reach_tol = 0.3;

% world
w_obs_min = 0.5 ; % minimum obstacle width [m] 
w_obs_max = 1.0 ; % maximum obstacle width [m]

n_dim = 2 ;
world_bounds = [0,15,0,8] ; % 2-D world
n_obs = 5; % number of obstacles

% agent
r_agent = 0.25; % m
v_max = 2.0; % m/s
w_max = 1.0; % rad/s

% planner
delta_v = 0.25;
delta_w = 0.25;

%% agent setup
x0 = zeros(4,1);

A = [1 0 dt 0; 
     0 1 0 dt;
     0 0 1 0;
     0 0 0 1];
 
B = [];
    
%% world setup
% start in bottom left and goal in upper right
start = [0.5;0.5]; goal = [world_bounds(2)-0.5;world_bounds(4)-0.5];

% generate uniformly-distributed random obstacles in middle portion of
% world
obs_bounds = world_bounds + [1.5,-1.5,0,0];
O_ctr = rand_in_bounds(obs_bounds,[],[],n_obs) ; % obstacle positions
O_wid = rand_range(w_obs_min,w_obs_max,[],[],2,n_obs) ; % obstacle widths

% obstacle zonotope array
obs = cell(1,n_obs) ;
for i = 1:n_obs
    % create generators from widths and randomly rotate 
    O_gen = O_wid(:,i) .* eye(n_dim) ;
    O_gen = rotation_matrix_2D(pi*rand()) * O_gen;
    obs{i} = ellipsotope(2,O_ctr(:,i),O_gen,[],[],{1,2}) ;
end

%% coloring parameters

color_obs = 'r' ; 
color_goal = 'g' ;
color_reach = 'b' ;
color_face_robot = [150 241 255] / 255 ; % light blue
color_edge_robot = [0 32 96] / 255 ; % dark blue

color_plan = [0.5 0.5 1] ; % blue
color_frs = [205 248 255] / 255 ; % light blue

%% plotting setup

figure(1) ; clf ; axis equal ; hold on ; grid on ;
axis(1.25.*(world_bounds))

% obstacles
for idx = 1:n_obs
    plot(obs{idx},'facecolor',color_obs,'edgecolor',color_obs,'facealpha',0.1);
end

% start positions
plot_path(start(1:n_dim),'go','markersize',15,'linewidth',1.5)

% goal areas
plot_path(goal,'gp','markersize',15,'markerfacecolor',color_goal)

%% RRT planning

num_RRT_iters = 1000;
iter_count = 0;

% params
RRT.step_dist = 0.5;
RRT.step_N = 20;
% tree variables
RRT.N_nodes = 0;
RRT.nodes = [];
RRT.nodes_parent = [];
RRT.parent_traj = []; % store trajectory from parent 

% initialize tree
disp('Initializing RRT Tree');
init_node = [start;0;0;0];
RRT.nodes = [RRT.nodes, init_node];
RRT.N_nodes = RRT.N_nodes + 1;

% grow tree
while iter_count < num_RRT_iters
    % pick random node from the tree
    idx_rand = rand_int(1,RRT.N_nodes);
    z_rand = RRT.nodes(:,idx_rand);
    
%     % sample random point
%     rand_pt = rand_in_bounds(world_bounds);
    % sample random control input
    u_rand = sample_control_input(z_rand,delta_v,delta_w,v_max,w_max);
    
%     % determine nearest node in RRT graph
%     node_distances = vecnorm(RRT.nodes - repmat(rand_pt,1,RRT.N_nodes));
%     [nearest_dist,nearest_idx] = min(node_distances) ;
%     nearest_node = RRT.nodes(:,nearest_idx) ;
%     
%     % create new node from nearest node in direction of sampled point
%     % (CURRENTLY IGNORES COLLISIONS)
%     diff = rand_pt - nearest_node;
%     heading = atan2(diff(2),diff(1));
%     new_node = nearest_node + RRT.step_dist * [cos(heading);sin(heading)];

    % compute trajectory from control input
    X = dubins_traj(z_rand(1:3),u_rand,dt,RRT.step_N);
    
    % check that trajectory does not intersect with obstacle or is out of bounds
    if ~is_in_obs(X,obs) && all(check_point_in_bounds(X(1:2,:),world_bounds))
        % add node to the tree
        new_node = [X(:,end);u_rand];
        RRT.nodes = [RRT.nodes, new_node];
        RRT.nodes_parent = [RRT.nodes_parent, idx_rand];
        RRT.N_nodes = RRT.N_nodes + 1;
        RRT.parent_traj = cat(3,RRT.parent_traj,X);
        
        % check if we reached the goal
        if norm(new_node(1:2) - goal) < goal_reach_tol
            disp('Goal Reached');
            break 
        end
        
        % update plots
        figure(1); hold on
        scatter(new_node(1),new_node(2));
        plot(X(1,:),X(2,:));
    end
    
    iter_count = iter_count + 1;
    if mod(iter_count,10) == 0
        disp(['iteration ',num2str(iter_count)])
    end
    
    % pause to allow plots to update
    pause(dt);
end

figure(1); hold on
% plot nodes
scatter(RRT.nodes(1,:),RRT.nodes(2,:))
% plot edges
for i = 2:length(RRT.nodes)
    traj = RRT.parent_traj(:,:,i-1);
    plot(traj(1,:),traj(2,:));
end

%% helper functions

% Given a trajectory and obstacle cell array (of ellipsotopes), determine 
% if any of the points of the trajectory lie in the obstacles
function check = is_in_obs(traj, obs)
    check = false;
    for i = 1:length(traj)
        point = traj(1:2,i);
        for j = 1:length(obs)
            if obs{j}.contains(point)
                check = true;
                return
            end
        end
    end
end

% Check if a trajectory lies within the world bounds
function check = in_bounds(traj, w_bounds)
    check = true;
    for i = 1:length(traj)
        point = traj(1:2,i);
        if ~check_point_in_bounds(point,w_bounds)
            check = false;
            return
        end
    end
end

% z - state
%  (x,y,theta,v,w)
function u = sample_control_input(z,delta_v,delta_w,v_max,w_max)
    % speed
    v_des_lo = max(z(4) - delta_v, 0);
    v_des_hi = min(z(4) + delta_v, v_max);
    v_des = rand_range(v_des_lo, v_des_hi, v_max, v_max/4);
    %v_des = v_max; % max speed for now
    % yaw rate
    w_des_lo = max(z(5) - delta_w, -w_max);
    w_des_hi = min(z(5) + delta_v, w_max);
    w_des = rand_range(w_des_lo,w_des_hi, z(5), w_max);
    % output
    u = [v_des;w_des];
end


function x_new = dubins_step(x,u,dt)
    x_new(1) = u(1)*cos(x(3))*dt + x(1);
    x_new(2) = u(1)*sin(x(3))*dt + x(2);
    x_new(3) = u(2)*dt + x(3);
end

function X_new = dubins_traj(x0,u,dt,N)
    X_new = zeros(3,N+1);
    x = x0; X_new(:,1) = x;
    for i = 2:N+1
        x = dubins_step(x,u,dt); 
        X_new(:,i) = x;
    end
end


% given the obstacles (as ellipsotopes) and world bounds, find a set of 
% initial locations in bounds which are sufficiently far from the obstacles 
% and from each other
function P_out = make_random_feasible_locations(n_agents,r_agents,O_buff,world_bounds)
    flag_P_feas = false ;
    
    while ~flag_P_feas
        % create random positions
        P_out = rand_in_bounds(world_bounds,[],[],n_agents) ;
        
        n_dim = size(P_out,1) ;
        
        % check position distances to obstacles
        obs_feas = true ;
        obs_buffer = 0.25 ;
        % for each position
        for i = 1:n_agents
            % for each obstacle
            for j = 1:length(O_buff)
                % convert position to zonotope with some buffering
                E_pos = ellipsotope(2,P_out(:,i),obs_buffer * eye(n_dim)) ;
                % check if they intersect
                if is_intersecting(E_pos,O_buff{j})
                    obs_feas = false ;
                    break
                end
            end
            if ~obs_feas
                break
            end
        end
        
        % check position distances to each other (greater than 5 agent
        % radiuses apart)
        D_to_self = dist_points_to_points(P_out,P_out) ;
        D_to_self(D_to_self == 0) = nan ; % ignore same points
        pos_feas = ~any(D_to_self(:) < 5*r_agents) ;
        
        % if positions are not too close to obstacle and not too close to
        % each other, break
        if obs_feas && pos_feas
            break
        end
    end
end

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

dt = 0.1;

buffer = 1;

w_obs_min = 0.5 ; % minimum obstacle width [m] 
w_obs_max = 1.0 ; % maximum obstacle width [m]

n_dim = 2 ;
world_bounds = 5.*[-1,1,-1,1] ; % 2-D world
n_obs = 10; % number of obstacles

r_agent = 0.25;
goal_reach_tol = 0.3;

%% agent setup
x0 = zeros(4,1);

A = [1 0 dt 0; 
     0 1 0 dt;
     0 0 1 0;
     0 0 0 1];
 
B = [];
    
%% world setup
% generate uniformly-distributed random obstacles
O_ctr = rand_in_bounds(world_bounds,[],[],n_obs) ; % obstacle positions
O_wid = rand_range(w_obs_min,w_obs_max,[],[],2,n_obs) ; % obstacle widths

% obstacle zonotope array
obs = cell(1,n_obs) ;
for i = 1:n_obs
    % create generators from widths and randomly rotate 
    O_gen = O_wid(:,i) .* eye(n_dim) ;
    O_gen = rotation_matrix_2D(pi*rand()) * O_gen;
    obs{i} = ellipsotope(2,O_ctr(:,i),O_gen,[],[],{1,2}) ;
end

% generate start and goal positions
locs = make_random_feasible_locations(2,r_agent,obs,world_bounds);
start = locs(:,1); goal = locs(:,2);

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

%% figure plot

% robot body 
G_robot = rotation_matrix_2D(pi*0.2) * diag([1 2]);
E_robot = ellipsotope(2,zeros(2,1),G_robot,[],[],{1,2});

% robot uncertainty
G_uncrt = rotation_matrix_2D(pi*-0.2) * diag([1 2]);
E_uncrt = ellipsotope(2,zeros(2,1),G_uncrt);

% mink sum
E_mink = E_robot + E_uncrt;

f = figure(2); 
subplot(1,2,1); axis equal; grid on;
plot(E_uncrt,'facecolor','r','edgecolor','r','facealpha',0.4)
xlabel('Uncertainty in $p\langle 1 \rangle$','Interpreter','latex'); ylabel('Uncertainty in $p\langle 2 \rangle$','Interpreter','latex');
lim = axis; axis(lim + [-1 1 -1 1]);

subplot(1,2,2); axis equal; grid on;
plot(E_robot,'facecolor','b','edgecolor','b','facealpha',1.0)
plot(E_mink,'facecolor','m','edgecolor','m','facealpha',0.4)
xlabel('$p\langle 1 \rangle$','Interpreter','latex'); ylabel('$p\langle 2 \rangle$','Interpreter','latex');
legend('Robot body','Sum of body and uncertainty');
lim = axis; axis(lim + [-1 1 -1 1]);

%save_figure_to_pdf(f,'../../.../figures/ellipsotopes/path_planning_uncertainty_sum.pdf')
save_figure_to_pdf(f,'path_planning_uncertainty_sum.pdf')
%% RRT planning

num_RRT_iters = 100;
iter_count = 0;

RRT.step_dist = 0.5;
RRT.N_nodes = 0;
RRT.nodes = [];
RRT.nodes_parent = [];

% initialize tree
disp('Initializing RRT Tree');
RRT.nodes = [RRT.nodes, start(1:2)];
RRT.N_nodes = RRT.N_nodes + 1;

% grow tree
while iter_count < num_RRT_iters
    % sample random point
    rand_pt = rand_in_bounds(world_bounds);
    
    % check if new point is in an obstacle
    if is_in_obs(rand_pt,obs)
        continue
    end
    
    % determine nearest node in RRT graph
    node_distances = vecnorm(RRT.nodes - repmat(rand_pt,1,RRT.N_nodes));
    [nearest_dist,nearest_idx] = min(node_distances) ;
    nearest_node = RRT.nodes(:,nearest_idx) ;
    
    % create new node from nearest node in direction of sampled point
    % (CURRENTLY IGNORES COLLISIONS)
    diff = rand_pt - nearest_node;
    heading = atan2(diff(2),diff(1));
    new_node = nearest_node + RRT.step_dist * [cos(heading);sin(heading)];
    
    % add node to the tree
    RRT.nodes = [RRT.nodes, new_node] ;
    RRT.nodes_parent = [RRT.nodes_parent, nearest_idx] ;
    RRT.N_nodes = RRT.N_nodes + 1 ;
    
    % check if we reached the goal
    if norm(new_node - goal) < goal_reach_tol
        disp('Goal Reached');
        break
    else
        iter_count = iter_count + 1;
        if mod(iter_count,10) == 0
            disp(['iteration ',num2str(iter_count)])
        end
    end
    
    % TODO: update plots
end

scatter(RRT.nodes(1,:),RRT.nodes(2,:))

%% helper functions

% Given a point and obstacle cell array (of ellipsotopes), determine if the
% point lies inside any of the obstacles
function check = is_in_obs(point, obs)
    check = false;
    for i = 1:length(obs)
        if obs{i}.contains(point)
            check = true;
            return
        end
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

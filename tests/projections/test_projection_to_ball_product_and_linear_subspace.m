%% description
% This script illustrates an iterative numerical method for finding points
% on the boundary of the intersection of a plane and a product of balls.
%
% See also: test_ball_slice_lemma.m,
% test_projection_to_constraint_and_ball.m
%
% NOTE: This doesn't yet work as of 5 Apr 2021 at 9:05 PM
%
% Authors: Shreyas Kousik
% Created: 5 Apr 2021
% Updated: 6 Ape 2021
clear;clc
%% user parameters
% rng seed
rng(0) ;

% p norm
p_norm = 2 ;

% index set
I = {1,[2,3]} ;

% linear constraint defining a plane a points for which Ax = b
A = [1 0 2] ;
b = 0.5 ;

% number of points to project
n_P = 1000 ;

% tolerance for getting points to constraint boundary
n_iter = 10 ;

%% automated from here
% get null space of linear constraint
t = pinv(A)*b ;
K = null(A) ;

% get dimensions of things
n_A = size(A,2) ; % dimension of ball space
n_K = size(K,2) ; % dimension of null space
n_I = length(I) ; % number of index subsets

if ~check_index_set_validity(I)
    error('Index set is not valid!')
end

%% get boundary of intersection of constraints
% sample a bunch of points in the space
P_orig = 2*rand(n_A,n_P) - 1 ;

% make random list of which norm to enforce per point
idxs_J_to_enf = rand_int(1,n_I,[],[],1,n_P) ;

% project points to linear constraint
F_proj = A*A'*pinv(A'*A) ;
P_proj = F_proj*(P_orig - t/2) ;
P_plane = P_orig - P_proj + repmat(t,1,n_P)/2 ;

%% iterate between projections for many points
% generate initial points
n_P = 1000 ;
P = 2*rand(n_A,n_P) - 1 ;

for idx = 1:n_iter
    % project points to ball
    P = project_points_to_ball_product(P,p_norm,I,idxs_J_to_enf) ;
    
    % project points to plane
    P = P - (F_proj*(P - t/2)) + t/2;
end

%% plotting setup
% create hyperplane
F_H = [1 2 3 4 1] ;
V_H = 2.*[K' ; -K'] + repmat(t',4,1) ;

% create ball product for patch
[F_E,V_E] = make_ball_product_for_patch(p_norm,I) ;

%% plotting
% create figure
fh = figure(1) ; clf ; axis equal ; hold on ; grid on ; view(3)

% plot ball product
h_E = patch('faces',F_E,'vertices',V_E,'facealpha',0.1','edgealpha',0,'facecolor','b') ;

% plot hyperplane
h_H = patch('faces',F_H,'vertices',V_H,'facealpha',0.1','edgealpha',0,'facecolor','r') ;

% plot final points
h_P = plot_path(P,'.','color','r','markersize',10) ;

legend([h_E,h_H,h_P],{'ball product','plane','final points'})

xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
set(gca,'fontsize',15)

save_figure_to_png(fh,'ellipsotope_proj_to_cons.png')
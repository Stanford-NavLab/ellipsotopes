%% description
% This script illustrates an iterative numerical method for finding points
% on the boundary of the intersection of a superellipsoid and a plane.
%
% See also: test_ball_slice_lemma.m
%
% Authors: Shreyas Kousik
% Created: 22 Mar 2021
% Updated: 5 Apr 2021
clear;clc
%% user parameters
% rng seed
rng(7) ;

% p-norm of ball
p = 6 ;

% linear constraint, Ax = b
A = [1 -1 2] ;
b = 0.5 ;

% number of points to project
n_P = 10 ;

% tolerance for getting points to constraint boundary
n_iter = 10 ;

%% automated from here
% get null space of linear constraint
t = pinv(A)*b ;
K = null(A) ;

% get dimensions of things
n_A = size(A,2) ; % dimension of ball space
n_K = size(K,2) ; % dimension of null space

%% get boundary of intersection of constraints
% sample a bunch of points in the space
P_orig = 2*rand(n_A,n_P) - 1 ;

% project points to linear constraint
F_proj = A*A'*pinv(A'*A) ;
P_proj = F_proj*(P_orig - t/2) ;
P_plane = P_orig - P_proj + repmat(t,1,n_P)/2 ;

%% plotting setup
% create hyperplane
F_H = [1 2 3 4 1] ;
V_H = 2.*[K' ; -K'] + repmat(t',4,1) ;

% create ball
[F_E,V_E] = make_superellipse_3D(p,1,zeros(3,1),500) ;

% plotting
% create figure
fh = figure(1) ; clf ; axis equal ; hold on ; grid on ; view(3)

% plot superellipsoid
h_E = patch('faces',F_E,'vertices',V_E,'facealpha',0.1','edgealpha',0,'facecolor','b') ;

% plot hyperplane
h_H = patch('faces',F_H,'vertices',V_H,'facealpha',0.1','edgealpha',0,'facecolor','r') ;

% plot points on plane
plot_path(P_plane,'k.','markersize',12)

%% iterate between projections for many points
% generate initial points
n_P = 1000 ;
P = 2*rand(n_A,n_P) - 1 ;
h_P_i = plot_path(P,'b.','markersize',1) ;

for idx = 1:n_iter
    % project points to ball
    N = vecnorm(P,p) ;
    P = P./repmat(N,n_A,1) ;
    
    % plot
    col = [idx/n_iter, 0, 1 - (idx/n_iter)] ;
    plot_path(P,'.','color',col,'markersize',2*idx)
    
    % project points to plane
    P = P - (F_proj*(P - t/2)) + t/2;
   
    % plot
    col = [idx/n_iter, 0, 1 - (idx/n_iter)] ;
    h_P = plot_path(P,'.','color',col,'markersize',2*idx) ;
    
    pause(0.1)
end

%% plotting cleanup
legend([h_E,h_H,h_P_i,h_P],...
    {'ellipsoid','plane',...
    'initial points','final points'})

xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
set(gca,'fontsize',15)

% save_figure_to_png(fh,'ellipsotope_proj_to_cons.png')
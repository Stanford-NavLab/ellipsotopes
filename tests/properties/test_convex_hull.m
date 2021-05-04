%% description
% This script tests computing the convex hull of two ellipsotopes.
%
% Authors: Adam Dai (and Shreyas sneakin in here too)
% Created: 19 Apr 2021
% Updated: 4 May 2021 (updated plotting to use ray tracing)
clear ; clc
%% user parameters
% rng seed
rng(0)

% ellipsotopes
p_norm = 2 ;

c_1 = zeros(2,1) ;
G_1 = eye(2) ;
A_1 = [] ;
b_1 = [] ;

c_2 = 2*rand(2,1) - 1 ;
G_2 = 2*rand(2) - 1 ;
A_2 = [] ;
b_2 = [] ;

%% automated from here
E1 = ellipsotope(p_norm,c_1,G_1,A_1,b_1);
E2 = ellipsotope(p_norm,c_2,G_2,A_2,b_2);

E_CH = convhull(E1,E2);

%% plotting
figure(1); clf; axis equal; hold on; grid on;

% plot etopes
% plot(E1);
% plot(E2);
% plot(E_CH);
plot_ray_tracing_2D(E1,100,'facecolor','r','edgecolor','r') ;
plot_ray_tracing_2D(E2,100,'facecolor','g','edgecolor','g') ;
plot_ray_tracing_2D(E_CH,200) ;
% plot_path(X(:,in_log),'r.')

%% testing if points are in/out of E_CH
% make grid of points
% n_X_1 = 50 ;
% n_X_2 = 25 ;
% X = make_grid_2D(3.*[-1,1,-1,1],n_X_1,n_X_2) ;
% n_X = size(X,2) ;
% in_log = false(1,n_X) ;
% 
% % for each point, test if it is inside the etope
% tic
% parfor idx_x = 1:n_X
%     x = X(:,idx_x) ;
%     chk = contains(E_CH,x) ;
%     
%     if chk
%         in_log(idx_x) = true ;
%     end    
% end
% toc
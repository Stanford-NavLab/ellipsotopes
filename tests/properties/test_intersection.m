%% description
% This script intersects 2 circular ellipsotopes to produce a constrained
% ellipsotope
%
% Authors: Adam Dai and Shreyas Kousik
% Created: shrug
% Updated: 15 Mar 2021
%
%% user parameters
% specify the two ellipsotopes
E1 = ellipsotope(2,[0;0],[1 0.5 -0.5; 0 0.866 0.866],[],[],{1,2,3});
E2 = ellipsotope(2,[0;0],diag([1;2]));

% or generate random ellipsotopes
gen_random_flag = true;

% whether or not to save flag
flag_save_figure = true;

%% automated from here
% perform the intersection 
E_int = E1 & E2;

%% plotting
h = figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot_ray_tracing_2D(E_int,200,'facecolor','m','edgecolor','m','facealpha',0.1);
plot_ray_tracing_2D(E1,200,'facecolor','r','edgecolor','r','facealpha',0.1);
plot_ray_tracing_2D(E2,200,'facecolor','b','edgecolor','b','facealpha',0.1);

set(gca,'fontsize',15)
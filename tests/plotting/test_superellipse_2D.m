%% description
% This script plots a 2-D superellipse defined by the equation
%
%   (x_1)^p + (x_2)^p = r^p
%
% In particular, we take a 2-D circle and dilate it by using the ratio
% of the 2-norm to the desired norm.
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2021
% Updated: nah
%
%% user parameters
% p-norm
p = 3 ;

% radius
r = 2 ;

% number of points to plot
n = 1000 ;

%% automated from here
% create points on unit circle
C = make_circle(r,n) ;

% get p-norm of points
P_C = vecnorm(C,p) ;
E = r.*C./P_C ;

% test function version
[F,V] = make_superellipse_2D(n,p,r,zeros(2,1)) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

patch('faces',F,'vertices',V,'facealpha',0.1) ;
plot_path(C,'b--')
plot_path(E,'b')

set_plot_linewidths(1.5)
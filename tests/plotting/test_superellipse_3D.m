%% description
% This script plots a 3-D superellipse defined by the equation
%
%   (x_1)^p + (x_2)^p + (x_3)^p = r^p
%
% In particular, we take a 3-D sphere and dilate it by using the ratio
% of the 2-norm to the desired norm.
%
% Authors: Shreyas Kousik
% Created: 05 Mar 2021
% Updated: 17 Mar 2021
%
%% user parameters
% p-norm
p = 4 ;

% radius
r = 0.6 ;

% approximate number of points to generate
n = 256 ;

%% automated from here
% adjust n
n_sphere = floor(sqrt(n-1)) ;

% create sphere
X = make_sphere(r,zeros(3,1),n) ;

% get vecnorms in the p-norm
N_X = vecnorm(X,p) ;

% dilate the sphere
S = r.*X./N_X ;

%% test functionized version
[F,V] = make_superellipse_3D(n,p,r,zeros(3,1)) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ; view(3)

plot_path(S,'b.')
plot_path(V','r.')
patch('faces',F,'vertices',V,'facealpha',0.1,'edgealpha',0.3)
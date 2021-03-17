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
% Updated: 15 Mar 2021
%
%% user parameters
% p-norm
p = 8 ;

% radius
r = 0.6 ;

% approximate number of points to generate
n = 1000 ;

%% automated from here
% adjust n
n_sphere = floor(sqrt(n-1)) ;

% create sphere
[X_1,X_2,X_3] = sphere(n_sphere) ;
X = [X_1(:), X_2(:), X_3(:)]' ;

% get vecnorms in the p-norm
N_X = vecnorm(X,p) ;

% dilate the sphere
S = r.*X./N_X ;

%% plotting setup
% get individual points
S_1 = S(1,:) ;
S_2 = S(2,:) ;
S_3 = S(3,:) ;

n_S = sqrt(length(S_1)) ;

% reshape points
S_1 = reshape(S_1,n_S,n_S) ;
S_2 = reshape(S_2,n_S,n_S) ;
S_3 = reshape(S_3,n_S,n_S) ;

% create patch
% [F,V] = surf2patch(S_1,S_2,S_3) ;

%% test functionized version
[F,V] = make_superellipse_3D(p,r,zeros(3,1),n) ;
figure(1) ; axis equal ; hold on ; grid on ; view(3) ;
patch('faces',F,'vertices',V,'facealpha',0.1)

%% plotting=
figure(1) ; clf ; axis equal ; hold on ; grid on ; view(3)

% plot_path(r.*X) ;
plot_path(S,'b.')
% plot_path(X_1,'r-')
% plot_path(X_2,'g-')
% plot_path(X_3,'k-')
patch('faces',F,'vertices',V,'facealpha',0.1,'edgealpha',0)
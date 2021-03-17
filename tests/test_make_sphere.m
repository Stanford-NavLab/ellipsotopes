%% description
% Make n points on a 3-D sphere of radius r centered at c \in \R^3. These
% points are evenly distributed over the surface of the sphere using
% Vogel's method as outlined in [1,2].
%
% References:
%   [1] https://bduvenhage.me/geometry/2019/07/31/generating-equidistant-vectors.html
%   [2] http://blog.marmakoide.org/?p=1
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2021
% Updated: 17 Mar 2021
%
%% user parameters
% radius
r = 2 ;

% center
c = rand(3,1) ;

% number of points
n = 200 ;

%% automated from here
% compute golden angle per http://blog.marmakoide.org/?p=1
ga = pi*(3 - sqrt(5)) ;

% create theta angles
T = ga * (0:(n-1)) ;

% create z values
Z = linspace(1 - (1/n), (1/n) - 1, n) ;

% create radii
R = sqrt(1 - Z.^2) ;

% create points at the given cylindrical coordinates
P = [R.*cos(T) ;
    R.*sin(T) ;
    Z] ;

% dilate points by r and shift to c
P = r.*P + repmat(c,1,n) ;

%% test functionized version
[F,V] = make_sphere(r,c,n) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ; view(3) ;

plot_path(P,'b.','markersize',12)
patch('faces',F,'vertices',V,'facealpha',0.1)
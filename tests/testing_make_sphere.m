%% description
%
%% user parameters
% radius
r = 2 ;

% center
c = zeros(2,1) ;

% number of points
n = 100 ;

%% automated from here
% get number of points on sphere
n_sph = n ; %floor(sqrt(n)) ;

% create uniform distribution of elevations
E = linspace(-pi/2,pi/2,n_sph) ;

% create distribution of azimuths that makes a roughly uniform distribution
% of points on the sphere; this should move "away" from 0 when E = -pi/2
% initially, "slow down" near E = 0, and then "speed up" again as E = pi/2
A = 2*pi*sin(sqrt(n_sph)*pi*E) ;

% create points at the given azimuths and elevations
P = [cos(A).*cos(E) ;
     -sin(A) ;
     -cos(A).*sin(E)] ;

%% study distribution of points
% distance from each point to each point
D = dist_points_to_points(P,P) ;
D = D(:) ;

% remove zeros
D(D==0) = [] ;
 
%% plotting
% plot sphere
figure(1) ; clf ; axis equal ; hold on ; grid on ; view(3) ;
plot_path(P,'b.')

% plot distribution of points on sphere
figure(2) ; clf ;
histogram(D)
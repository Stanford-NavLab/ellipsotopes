%% description
% This script visualizes our lemma that the intersection of a plane with a
% 2-norm ball in n dimensions is the affine map of a 2-norm ball in n-1
% dimensions.
%
% Authors: Shreyas Kousik
% Created: shrug
% Updated: 22 Mar 2021
%
%% user parameters
% hyperplane in 3-D
A = [1 1 1] ;
b = 0.5 ;

%% automated from here
% get affine map to hyperplane
t = pinv(A)*b ; % center of intersection area
c = sqrt(1 - vecnorm(t)^2) ; % radius of intersected ball

K = null(A) ;
T = c*K ;

%% plotting setup
% create unit sphere
[F_E,V_E] = make_ellipsoid_for_patch(1,1,1) ;

% create hyperplane
F_H = [1 2 3 4 1] ;
V_H = 2.*[K' ; -K'] + repmat(t',4,1) ;

% create lower-dimensional ball
n = 100 ;
[F_C,V_C] = make_circle(1,n) ;
V_C = (T*V_C' + repmat(t,1,n))' ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ; view(3)

% plot ellipsoid
patch('faces',F_E,'vertices',V_E,'facealpha',0.1','edgealpha',0,'facecolor','b')

% plot hyperplane
patch('faces',F_H,'vertices',V_H,'facealpha',0.1','edgealpha',0,'facecolor','r')

% plot intersection
patch('faces',F_C,'vertices',V_C,'facealpha',0.1','edgealpha',1,'facecolor','g')

% plot arrow to t
plot_arrow(t)
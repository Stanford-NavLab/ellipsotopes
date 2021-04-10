%% description
% This script tests the containment lemma where least-squares is used to
% determine if a point is in/out of an ellipsotope
%
% Authors: Shreyas Kousik
% Created: 5 Mar 2021
% Updated: nah
%
%% user parameters
% ellipsotope params
c = [0;0] ;
G = [1 2 ; -1 1] ;
p_norm = 6 ;

% constrained ellipsotope


%% automated from here
% make ellipsotope
E = ellipsotope(p_norm,c,G) ;

% make grid on points
P = make_grid_2D(3*[-1,1,-1,1],50,50) ;

% test if point is inside the ellipsotope
A = E.generators ;
b = P - E.center ;
t = pinv(A)*b ;

flag_in = all(t <= 1,1) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E)
plot_path(P(:,flag_in),'b.')
plot_path(P(:,~flag_in),'r.')
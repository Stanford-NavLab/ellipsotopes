%% description
% This script is a sanity check to make sure I remember linear algebra...
%
% Authors: Shreyas Kousik
% Created: 9 Apr 2021
% Updated: not yet
clear ; clc ;
%% user parameters
A = [1 1 1 ;
    -1 -1 1] ;
b = 0.1*ones(2,1) ;

n_P = 100 ;

%% automated from here
% get dimension
n_con = size(A,2) ;

% make points
P = 2*rand(n_con,n_P) - 1;

% get null space and a feasible point
K = null(A) ;
t = A\b ;

% project points onto null space
P_dots = K'*(P - repmat(t,1,n_P)) ;
P_proj = K*P_dots ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ; view(3)

% plot linear subspace planes
plot_planes_3D(A,b)

% plot original points
plot_path(P,'r.','markersize',8)

% plot projected ooints
plot_path(P_proj,'b.','markersize',8)
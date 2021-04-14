%% description
% This script tests generating line segments spanning an ellipsotope by
% intersecting it with two halfplanes facing opposite directions. Right now
% this is just a wild idea, but it might really speed up plotting these
% objects!
%
% Authors: Shreyas Kousik
% Created: 13 Apr 2021
% Updated: nuuu
clear ; clc ;
%% user parameters
% rng seed
rng(0)

% ellipsotope
p_norm = 4 ;
c = rand(2,1) ;
G = rand(2) ;
A = [] ;
b = [] ;
I = {[1,2]} ;

% number of directions to generate (DON'T MAKE THIS TOO LARGE)
n_U = 5 ;

%% automated from here
% generate a bunch of directions
U = [rotation_matrix_2D(pi/2)*G, 2*rand(2,n_U) - 1] ;

figure(1) ; clf ; axis equal ; grid on ; hold on ;

for idx = 1:n_U
    % get direction
    u = U(:,idx)' ;
    
    % perform intersections
    [c_temp,G_temp,A_temp,b_temp,I_temp] = intersect_ellipsotope_with_halfspace(zeros(2,1),G,A,b,I,u,0) ;
    [c_int,G_int,A_int,b_int,I_int] = intersect_ellipsotope_with_halfspace(c_temp,G_temp,A_temp,b_temp,I_temp,-u,0) ;

    % get the nullspace of the constraint
    K = null(A_int) ;
    t = A_int\b_int ;
    
    % plot intersection with pair of halfplanes
    plot_ellipsotope_utility(p_norm,c,G_int,A_int,b_int,I_int,'edgecolor','r')
    
end
%% plotting
% figure(1) ; clf ; axis equal ; grid on ; hold on ;

% plot original ellipsotope
plot_ellipsotope_utility(p_norm,c,G,A,b,I)

% % plot intersection with pair of halfplanes
% plot_ellipsotope_utility(p_norm,c,G_int,A_int,b_int,I_int,'edgecolor','r')

%% helper functions
function [c_int,G_int,A_int,b_int,I_int] = intersect_ellipsotope_with_halfspace(c,G,A,b,I,h,f)
% get sizes of things
n_gen = size(G,2) ;
n_con = size(A,1) ;

    d = 0.5*(sum(abs(h*G),2) + f - h*c) ;
    c_int = c ;
    G_int = [G, zeros(2,1)] ;
    A_int = [A, zeros(n_con,1) ; h*G, d] ;
    b_int = [b ; f - h*c - d] ;
    I_int = [I, {n_gen + 1}] ;
end
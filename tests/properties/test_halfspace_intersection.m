%% description
% This script tests intersecting an ellipsotope with a halfplane -- AND IT
% WORKS!
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
c = zeros(2,1) ;
G = eye(2) ;
A = [] ;
b = [] ;
I = {1:2} ;

% halfplane {x | h*x <= f}
h = [1 0] ;
f = 0.3 ;

%% automated from here
% get sizes of things
n_gen = size(G,2) ;
n_con = size(A,1) ;

% perform intersection
d = 0.5*(sum(abs(h*G),2) + f - h*c) ;
G_S = [G, zeros(2,1)] ;
A_S = [A, zeros(n_con,1) ; h*G, d] ;
b_S = [b ; f - h*c - d] ;
I_S = [I, {n_gen + 1}] ;

%% plotting
figure(1) ; clf ; axis equal ; grid on ; hold on ;

plot_ellipsotope_utility(p_norm,c,G,A,b,I)
plot_ellipsotope_utility(p_norm,c,G_S,A_S,b_S,I_S)
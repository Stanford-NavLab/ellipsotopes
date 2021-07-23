%% description
% This script tests sampling points from an ellipsotope
%
% Authors: Shreyas Kousik
% Created: 23 July 2021
clear ; clc ;
%% user parameters
% rng seed
rng(0)

% ellipsotope properties
p_norm = 2 ;
n_dim = 2 ;
n_gen = 6 ;
n_con = 2 ;
n_I = 5 ;

% number of points to sample
n_P = 100 ;

%% automated from here
% make the etope
E = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con,n_I) ;

P = sample_from_ellipsotope(E,n_P) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E)
plot_path(P,'r.')

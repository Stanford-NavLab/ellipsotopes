%% description
% This script tests Proposition 4 of [1] for ellipsotopes; turns out, it
% doesn't work! Ha ha ha ha
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% Authors: Shreyas Kousik
% Created: 16 Mar 2022
% Updated: --
clear ; clc ;
%% user parameters
% rng seed

% etope specs
c = zeros(2,1) ;
G = [1 0 1 ; 1 2 -1] ;
A = [-1 1 -1] ;
b = 2 ;
I = {[1],[2,3]} ;

% n_dim = 2 ;
% n_gen = 13 ;
% n_con = 3 ;

%% automated from here
% E = make_random_ellipsotope(2,n_dim,n_gen,n_con) ;
E = ellipsotope(2,c,G,A,b,I) ;

% get interval containing intersection of {Ax=b} and infty-norm ball
% product space (i.e. hypercube)
[int_E,int_R,int_R_rescaled,xi_m,xi_r] = estimate_interval_ball_product_intersection(E) ;

% rescale constraints of E using (24) in [1], Prop. 4
[p,c,G,A,b,I] = get_properties(E) ;
G_rsc = G*diag(xi_r) ;
c_rsc = c + G*xi_m ;
A_rsc = A*diag(xi_r) ;
b_rsc = b - A*xi_m ;

% make new tope
E_rsc = ellipsotope(2,c_rsc,G_rsc,A_rsc,b_rsc,I) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(E)
plot(E_rsc,'color','r','linestyle','--','linewidth',3)

make_plot_pretty()
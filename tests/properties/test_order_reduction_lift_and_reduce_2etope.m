%% description
% This script tests the lift-and-reduce strategy from the constrained
% zonotope paper [1]. In particular we leverage the fact that
% 2-ellipsotopes can be order-reduced in a nice way.
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% Authors: Shreyas Kousik
% Created: 3 Mar 2022
% Updated: 10 Mar 2022 (functionized the lift-and-reduce strategy)

%% user parameters
% rng seed
rng(0) ;

% original tope
n_dim = 2 ;
n_gen = 30 ;
n_con = 3 ;
n_I = 3 ;

%% automated from here
% make original etope
[E,c,G,A,b,I] = make_random_ellipsotope(2,n_dim,n_gen,n_con,n_I) ;

E_rdc = reduce_2_etope_to_minimal_exact_rep(E) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(E)
plot(E_rdc,'color','r','linestyle','--','linewidth',3)

make_plot_pretty()
%% description
% Tests out the identify_component_zonotope function
%
% Authors: Shreyas Kousik
% Created: 14 Mar 2022
% Updated: 15 Mar 2022

clear ; clc ;
%% user parameters
% rng seed
rng(1)

% etope specs
n_dim = 2 ;
n_gen = 15 ;
n_con = 3 ;
n_I = 8 ;

%% automated from here
E = make_random_ellipsotope(2,n_dim,n_gen,n_con,n_I) ;
E = reduce_2_etope_to_minimal_exact_rep(E) ;

%% identify comp zono
[E_reorg,idx_start_Z_I,idx_start_Z] = identify_component_zonotope(E) ;

%% plotting
% figure(1) ; clf ; axis equal ; hold on ; grid on ;
% plot(E)
% plot(E_reorg,'color','r','linestyle','--')
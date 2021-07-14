%% description
% This script demonstrates reducing an ellipsotope's number of generators
% by identifying component ellipsotopes.
%
% Authors: Shreyas Kousik
% Created: 13 July 2021
% Updated: not yet
clear ; clc ;

%% automated from here
% make ellipsotope
rng(0)
E_comp_1 = make_random_ellipsotope(2,2,4,1,1) ;
E_comp_2 = make_random_ellipsotope(2,2,4,1,1) ;
E_other = make_random_ellipsotope(2,2,5,2,2) ;
E = E_comp_1 + E_other + E_comp_2 ;

% get propz
[p,c,G,A,b,I,n_dim,n_gen,n_con,n_I] = get_properties(E) ;

% create copies of propz
G_old = G ;
A_old = A ;
b_old = b ;
I_old = I ;

% identify component ellipsotopes
[idxs,log_idxs,E_reorg,E_other,E_comp_cell] = E.identify_component_ellipsotopes() ;

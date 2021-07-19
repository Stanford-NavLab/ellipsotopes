%% description
% This script demonstrates reducing an ellipsotope's number of generators
% by identifying component ellipsotopes.
%
% Authors: Shreyas Kousik
% Created: 13 July 2021
% Updated: 14 July 2021
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

% reduce each component ellipsotope (since we know it's constrained = basic!)
n_comp = length(E_comp_cell) ;
for idx = 1:n_comp
    E_comp_cell{idx} = reduce_constrained_2_etope(E_comp_cell{idx}) ;
end

%% reassemble original etope from component+other topes
E_reassembled = E_other ;
for idx = 1:length(E_comp_cell)
    E_reassembled = E_reassembled + E_comp_cell{idx} ;
end

%% use MVOE order reduction for component etopes
E_cell_out = reduce_component_2_etopes(E_comp_cell) ;

E_rdc = E_cell_out{1} ;

for idx = 2:length(E_cell_out)
    E_rdc = E_rdc + E_cell_out{idx} ;
end

E_rdc = E_rdc + E_other ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E)
plot(E_reassembled,'color','r','linewidth',3,'linestyle','--')
plot(E_rdc,'color','g')

legend('original','reassembled','reduced')

set_plot_fontsize(15)


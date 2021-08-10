%% description
% This script identifies component ellipsotopes by checking for block
% diagonal entries in the constraint matrix according to the index set.
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


% identify component ellipsotopes
[idxs,log_idxs,E_reorg,E_other,E_comp_cell] = identify_component_ellipsotopes(E) ;

%% reassemble original etope from component+other topes
E_reassembled = E_other ;
for idx = 1:length(E_comp_cell)
    E_reassembled = E_reassembled + E_comp_cell{idx} ;
end

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E)
plot(E_reassembled,'color','r','linewidth',3,'linestyle','--')

legend('original','reassembled')

set_plot_fontsize(15)


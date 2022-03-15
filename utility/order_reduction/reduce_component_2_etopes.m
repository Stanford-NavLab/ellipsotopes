function E_cell_out = reduce_component_2_etopes(E_cell_in,n_rdc)
% E_cell_out = reduce_component_2_etopes(E_cell,n_rdc)
%
% Given a list of basic 2-ellipsotopes E_cell, create a minimum-volume
% outer ellipsoid to eliminate n_rdc of them. By default, n_rdc is 1.
%
% The output is a list of basic 2-ellipsotopes that has n_rdc fewer etopes
% than the original input list.
%
% Authors: Shreyas Kousik
% Created: 14 July 2021
% Updated: 14 Mar 2022 (fixed how n_rdc is used)

    % set default number of topes to reduce
    if nargin < 2
        n_rdc = 1 ;
    end
    
    n_gen_in = get_total_number_of_generators(E_cell_in) ;
    n_gen_des = n_gen_in - n_rdc ; % desired number of generateurs

    E_cell_out = E_cell_in ;
    n_gen = n_gen_in ;

    % reduce reduce reduce!
    if ~isempty(E_cell_out)
        while n_gen > n_gen_des
            E_cell_out = reduce_helper_function(E_cell_out) ;
            n_gen = get_total_number_of_generators(E_cell_out) ;
        end
    end
end

function E_cell_out = reduce_helper_function(E_cell)
% E_cell_out = reduce_helper_function(E_cell_in)
%
% This helper function reduces the length of E_cell_in by 1. For now, this
% is a super inefficient implementation, but it'll do, pig, it'll do.
%
% Authors: Shreyas Kousik
% Created: 14 July 2021
% Updated: 13 Mar 2022 (updated to use new heuristic)

    % number of etopes
    n_E = length(E_cell) ;
    
    % create all (i,j) pairs
    combs = combinator(n_E,2,'c') ;
    n_idx = size(combs,1) ;

    % iterate over pairs to compute statistic
    vols_heur = nan(n_E,1) ;
    for idx = 1:n_idx
        % get ellipsoids to compare
        idx_ij = combs(idx,:) ;
        G_i = E_cell{idx_ij(1)}.generators ;
        G_j = E_cell{idx_ij(2)}.generators ;

        % compute heuristic
        Q_i = inv(pinv(G_i)'*pinv(G_i)) ;
        Q_j = inv(pinv(G_j)'*pinv(G_j)) ;
        Q_MVOE = 2*(Q_i + Q_j) ; % assume bt = 1 ;
        vols_heur(idx) = 1/det(inv(Q_MVOE)) ;
    end

    % get the pair of topes that minimize the heuristic
    [~,idx_min] = min(vols_heur) ;
    idx_ij = combs(idx_min,:) ;
    idx_i = idx_ij(1) ;
    idx_j = idx_ij(2) ;

    % extract the two topes
    E_i = E_cell{idx_i} ;
    E_j = E_cell{idx_j} ;
    E_cell(idx_ij) = []  ; % delete the topes from the input

    % compute new MVOE etope
    c_MVOE = E_i.center + E_j.center ;
    G_MVOE = make_MVOE_generator_matrix(E_i.generators,E_j.generators) ;
    E_MVOE = ellipsotope(2,c_MVOE,G_MVOE) ;
    
    % create output
    E_cell_out = [E_cell, {E_MVOE}] ;
end

function n_gen = get_total_number_of_generators(E_cell)
% n_gen = get_total_number_of_generators(E_cell)
%
%
    n_gen = 0 ;
    for idx = 1:length(E_cell)
        n_gen = n_gen + E_cell{idx}.n_generators ;
    end
end
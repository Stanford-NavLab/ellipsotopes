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
% Updated: nup

% set default number of topes to reduce
if nargin < 2
    n_rdc = 1 ;
end

E_cell_out = E_cell_in ;

for idx = 1:n_rdc
    E_cell_out = reduce_helper_function(E_cell_out) ;
end

end

function E_cell_out = reduce_helper_function(E_cell_in)
% E_cell_out = reduce_helper_function(E_cell_in)
%
% This helper function reduces the length of E_cell_in by 1. For now, this
% is a super inefficient implementation, but it'll do, pig, it'll do.
%
% Authors: Shreyas Kousik
% Created: 14 July 2021
% Updated: nup

    % number of etopes
    n_E = length(E_cell_in) ;

    % set up to save generator matrices
    G_cell = cell(1,n_E) ;

    % compute ellipsoid axes for each etope
    V = cell(1,n_E) ;
    for idx = 1:n_E
        % get current etope
        E_idx = E_cell_in{idx} ;
        G_idx = E_idx.generators ;
        G_cell{idx} = G_idx ;

        G_inv = pinv(G_idx) ;
        Q = G_inv'*G_inv ;
        [V_idx,E_idx] = eig(inv(Q)) ;

        % get axis for largest eigenvalue
        E_idx = diag(E_idx) ;
        [~,idx_max] = max(E_idx) ;

        V{idx} = V_idx(:,idx_max) ;
    end

    % create all (i,j) pairs
    combs = combinator(n_E,2,'c') ;
    n_idx = size(combs,1) ;

    % iterate over pairs to compute statistic
    rdc_heur = nan(n_E,1) ;
    for idx = 1:n_idx
        % get ellipsoids to compare
        idx_ij = combs(idx,:) ;

        % get longest axes
        v_i = V{idx_ij(1)} ;
        v_j = V{idx_ij(2)} ;

        % compute heuristic
        rdc_heur(idx) = abs(v_i'*v_j) ;
    end

    % get the pair of topes that max the heuristic
    [~,idx_max] = max(rdc_heur) ;
    idx_ij = combs(idx_max,:) ;
    idx_i = idx_ij(1) ;
    idx_j = idx_ij(2) ;

    % extract the two topes
    E_i = E_cell_in{idx_i} ;
    E_j = E_cell_in{idx_j} ;
    E_cell_in(idx_ij) = []  ; % delete the topes from the input

    % compute MVOE generator matrix
    c_MVOE = E_i.center + E_j.center ;
    G_MVOE = make_MVOE_generator_matrix(G_cell{idx_i},G_cell{idx_j}) ;
    E_MVOE = ellipsotope(2,c_MVOE,G_MVOE) ;
    
    % create output
    E_cell_out = [E_cell_in, {E_MVOE}] ;
end
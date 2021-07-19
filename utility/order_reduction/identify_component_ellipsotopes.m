function [idxs,log_idxs,E_reorg,E_other,E_comp_cell] = identify_component_ellipsotopes(E,flag_reorg)
% idxs = identify_component_ellipsotopes(E)
% idxs = identify_component_ellipsotopes(E,reorganize_topes_flag)
%
% Find all component ellipsotopes a given ellipsotope. The output indices
% correspond to the etope's index set. By default, this method also
% reorganizes the ellipsotope's data so that the component ellipsotopes
% are the last few index subsets.
%
% Authors: Shreyas Kousik
% Created: 13 July 2021
% Updated: 19 July 2021

    % get properties
    [~,~,~,A,~,I] = get_properties(E) ;
    n_I = length(I) ;

    % initialize second output
    log_idxs = true(1,n_I) ;

    if nargin < 2
        flag_reorg = true ;
    end

    %% identify component etopes
    % if constraint set is empty then every index subset is a component etope
    if isempty(A)
        idxs = 1:n_I ;
    else
        % initialize output
        idxs = [] ;

        % iterate through index set and check for component etopes
        for idx = 1:n_I
            % get current indices
            J = I{idx} ;

            % get all constraints with nonzero values at the current indices
            A_row_log = any(A(:,J) ~= 0,2) ;

            % check if all other entires in the current constraints are zeros
            A_temp = A ;
            A_temp(:,J) = [] ;
            if all(all(A_temp(A_row_log,:) == 0,2))
                idxs = [idxs, idx] ;
            else
                log_idxs(idx) = false ;
            end
        end
    end

    %% reorganize so the component etopes are the last ones
    if flag_reorg
        % get propz
        [p,c,G,A,b,~,n_dim] = get_properties(E) ;

        % we'll make a new generator matrix, constraints, and index set for
        % the component etopes and the other ones
        G_other = [] ;
        A_other = [] ;
        I_other = [] ;
        n_other = 0 ;

        G_comp = [] ;
        A_comp = [] ;
        I_comp = [] ;
        n_comp = 0 ;
        
        E_comp_cell = [] ; 

        for idx = 1:n_I
            J = I{idx} ;

            G_idx = G(:,J) ;
            A_idx = A(:,J) ;
            n_G_idx = size(G_idx,2) ;


            if log_idxs(idx)
                % component etope!
                I_temp = {1:n_G_idx} ;
                I_idx = shift_index_set(I_temp,n_comp) ;
                n_comp = n_comp + n_G_idx ;

                G_comp = [G_comp, G_idx] ;
                A_comp = [A_comp, A_idx] ;
                I_comp = [I_comp, I_idx] ;
                
                E_comp_idx = ellipsotope(p,zeros(n_dim,1),G_idx,A_idx,b,I_temp) ;
                E_comp_cell = [E_comp_cell, {E_comp_idx}] ;                
                
            else
                % other etope!
                I_idx = shift_index_set({1:n_G_idx},n_other) ;
                n_other = n_other + n_G_idx ;

                G_other = [G_other, G_idx] ;
                A_other = [A_other, A_idx] ;
                I_other = [I_other, I_idx] ;
            end
        end

        % after iterating through all the index subsets, shift the
        % component etopes to the "end" of the index set
        if ~isempty(I_comp)
            I_comp = shift_index_set(I_comp,n_other) ;
        end

        % reassemble the etope
        G_new = [G_other, G_comp] ;
        A_new = [A_other, A_comp] ;
        I_new = [I_other, I_comp] ;

        E_reorg = ellipsotope(p,c,G_new,A_new,b,I_new) ;
        E_other = ellipsotope(p,c,G_other,A_other,b,I_other) ;
    end
end
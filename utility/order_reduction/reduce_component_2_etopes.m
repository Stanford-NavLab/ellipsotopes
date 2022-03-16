function E = reduce_component_2_etopes(E,n_rdc)
% E_rdc = reduce_component_2_etopes(E,n_gen_to_remove)
%
% Given a 2-ellipsotope E, try to reduce n_rdc generators by identifying
% and MVOEing component ellipsoids.
%
% See also: reduce_component_2_etopes_cell.m,
%           test_order_reduction_component_ellipsotopes.m
%
% Authors: Shreyas Kousik
% Created: 14 July 2021
% Updated: 14 Mar 2022 (fixed how n_rdc is used)

    %% setup
    % get properties
    if E.p_norm ~= 2
        error('This only works for 2-ellipsotopes!')
    end
    
    n_dim = E.dimension ;
    
    if nargin < 2
        n_rdc = 1 ;
    end

    %% identify component ellipsotopes
    [~,~,~,E_other,E_comp_cell] = identify_component_ellipsotopes(E) ;

    %% if there are any component etopes, reduce them
    if ~isempty(E_comp_cell)
        % reduce each component ellipsotope to an ellipsoid using the
        % fact that constrained = basic for ellipsotopes
        n_comp = length(E_comp_cell) ;
        for idx = 1:n_comp
            E_comp_cell{idx} = reduce_constrained_2_etope(E_comp_cell{idx}) ;
        end

        % use MVOE order reduction for component etopes
        if length(E_comp_cell) > 1
            E_cell_out = reduce_component_2_etopes_cell(E_comp_cell,n_rdc) ;
        else
            E_cell_out = E_comp_cell ;
        end

        E = E_other ;

        if ~isempty(E_cell_out)
            for idx = 1:length(E_cell_out)
                E = E + E_cell_out{idx} ;
            end
        end
    end
end
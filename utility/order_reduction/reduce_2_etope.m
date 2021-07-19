function E_rdc = reduce_2_etope(E,n_rdc)
% E_rdc = reduce_2_etope(E)
%
% This function implements the 2-ellipsotope order reduction strategy in
% the ellipsotope paper.
%
% Authors: Shreyas Kousik
% Created: 19 July 2021
% Updated: not yet

    % sanity check
    if E.p_norm ~= 2
        error('This function only works for 2-ellipsotopes!')
    end
    
    % set default number of topes to reduce
    if nargin < 2
        n_rdc = 1 ;
    end
    
    % make sure n_rdc wouldn't make the ellipsotope smaller than the system
    % dimension
    if (E.n_generators - n_rdc) < E.n_dimension
        n_rdc = E.n_generators - E.n_dimension ;
    end

    % identify component ellipsotopes
    [~,~,~,E_other,E_comp_cell] = identify_component_ellipsotopes(E) ;

    % reduce each component ellipsotope (since we know it's constrained = basic!)
    n_comp = length(E_comp_cell) ;
    for idx = 1:n_comp
        E_comp_cell{idx} = reduce_constrained_2_etope(E_comp_cell{idx}) ;
    end

    % use MVOE order reduction for component etopes
    E_cell_out = reduce_component_2_etopes(E_comp_cell,n_rdc) ;

    E_rdc = E_other ;

    if ~isempty(E_cell_out)
        for idx = 1:length(E_cell_out)
            E_rdc = E_rdc + E_cell_out{idx} ;
        end
    end
end
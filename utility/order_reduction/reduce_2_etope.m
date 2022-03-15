function E = reduce_2_etope(E,n_rdc)
% E_rdc = reduce_2_etope(E)
%
% This function implements the 2-ellipsotope order reduction strategy in
% the ellipsotope paper. If the second input, n_rdc, is specified, then we
% try to get rid of that many generators.
%
% Authors: Shreyas Kousik
% Created: in days of yore
% Updated: 14 Mar 2022 (leveraging new proposition)

%% setup
    % sanity check
    if E.p_norm ~= 2
        error('This function only works for 2-ellipsotopes!')
    end
    
    % set default number of topes to reduce
    if nargin < 2
        n_rdc = 1 ;
    end
    
    % desired number of gennies
    n_gen_orig = E.n_generators ;
    n_des = n_gen_orig - n_rdc ;
    
%% get into minimal form
    E = reduce_2_etope_to_minimal_exact_rep(E) ;
    n_gen = E.n_generators ;
    
%% try using component ellipsotopes
    if n_gen > (n_gen_orig - n_rdc)
        % identify component ellipsotopes
        [~,~,~,E_other,E_comp_cell] = identify_component_ellipsotopes(E) ;
        
        if ~isempty(E_comp_cell)
            % reduce each component ellipsotope to an ellipsoid using the
            % fact that constrained = basic for ellipsotopes
            n_comp = length(E_comp_cell) ;
            for idx = 1:n_comp
                E_comp_cell{idx} = reduce_constrained_2_etope(E_comp_cell{idx}) ;
            end

            % use MVOE order reduction for component etopes
            if length(E_comp_cell) > 1
                E_cell_out = reduce_component_2_etopes(E_comp_cell,n_rdc) ;
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
        
        n_gen = E.n_generators ;
    end
    
%% try removing constraints so the minimal rep thing will get us low enough
    n_con = E.n_constraints ;
    while (n_gen > n_des) && (n_con > 0)
        E = reduce_etope_constraint_and_generator(E) ;
        E = reduce_2_etope_to_minimal_exact_rep(E) ;
        n_gen = E.n_generators ;
        n_con = E.n_constraints ;
    end

% %% try lift-and-reduce
%     if n_gen > n_des
%         E = reduce_2_etope_lift_and_reduce(E,n_gen - n_des) ;
%     end
%     
%     n_gen = E.n_generators ;
%     
% %% try component zonotopes?

    %% sanity check
    n_gen = E.n_generators ;
    if n_gen > n_des
        warning(['Still gotta write how to reduce component zonotopes! ',...
            'So the output of this function might not be as reduced as ',...
            'you desire.'])
    end
end
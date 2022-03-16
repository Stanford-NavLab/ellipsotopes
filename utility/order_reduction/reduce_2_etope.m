function E = reduce_2_etope(E,n_gen_to_remove,n_con_to_keep)
% E_rdc = reduce_2_etope(E)
% E_rdc = reduce_2_etope(E,n_gen_to_remove,n_con_to_keep)
%
% Inputs:
%   E - ellipsotope
%
%   n_gen_to_remove - number of generators to (attempt to) remove
%
%   n_con_to_keep - minimum number of constraints to keep
%
% Outputs:
%   E_rdc - reduced ellipsotope
%
% This function implements the 2-ellipsotope order reduction strategy in
% the ellipsotope paper. If the second input, n_rdc, is specified, then we
% try to get rid of that many generators.
%
% Authors: Shreyas Kousik
% Created: in days of yore
% Updated: 15 Mar 2022

%% setup
    % sanity check
    if E.p_norm ~= 2
        error('This function only works for 2-ellipsotopes!')
    end
    
    % set default number of topes to reduce
    if nargin < 2
        n_gen_to_remove = 1 ;
    end
    
    if nargin < 3
        n_con_to_keep = 3 ; % threeee is a magic number
    end
    
    if nargin < 4
        flag_force_reduce = true ;
    end
    
    % desired number of gennies
    n_gen_orig = E.n_generators ;
    n_des = n_gen_orig - n_gen_to_remove ;
    
%% get into minimal form
    E = reduce_2_etope_to_minimal_exact_rep(E) ;
    n_gen = E.n_generators ;
    
%% try using component ellipsotopes
    if n_gen > (n_gen_orig - n_gen_to_remove)
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
                E_cell_out = reduce_component_2_etopes(E_comp_cell,n_gen_to_remove) ;
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
    while (n_gen > n_des) && (n_con > n_con_to_keep)
        E = reduce_etope_constraint_and_generator(E) ;
        E = reduce_2_etope_to_minimal_exact_rep(E) ;
        n_gen = E.n_generators ;
        n_con = E.n_constraints ;
    end

%% try lift-and-reduce
    if n_gen > n_des
        E = reduce_2_etope_lift_and_reduce(E,n_gen - n_des) ;
    end
    
    n_gen = E.n_generators ;
    
%% try component (constrained) zonotopes
    if n_gen > n_des
        % lift!
        E_l = lift(E) ;
        [E_l,idx_Z_I,idx_Z] = identify_component_zonotope(E_l) ;
        [p,c_l,G_l,~,~,I_l] = get_properties(E_l) ;

        % figure out how much to reduce the zonotope generators
        n_rdc = n_gen - n_des ;
        
        % get the lifted generator matrix of the component zonotope
        G_rdc = G_l(:,idx_Z:n_gen) ;
        [n_r,n_c] = size(G_rdc) ;
        
        % if G_rdc is "wide" then it can be reduced, otherwise we need to
        % pop some generators...
        
        if n_c > n_r
            G_rdc = reduce_zonotope_Chischi(G_rdc,n_rdc) ;

            % create a new index set for the reduced generator matrix
            I_rdc = I_l(1:(idx_Z_I-1)) ;
            I_rdc = [I_rdc, num2cell((1:size(G_rdc,2))+idx_Z)] ;

            % reconstruct tope
            E_l = ellipsotope(p,c_l,G_rdc,[],[],I_rdc) ;
            E = drop(E_l,n_dim) ;
        end
    end

%% final check
    n_gen = E.n_generators ;
    if n_gen > n_des
        warning(['Still gotta write how to reduce as much as desired! ',...
            'So the output of this function might not be as reduced as ',...
            'you wanted :('])
    end
end
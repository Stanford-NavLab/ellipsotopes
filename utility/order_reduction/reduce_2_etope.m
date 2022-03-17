function E = reduce_2_etope(E,n_gen_to_remove,n_con_to_keep,flag_force_reduce)
% E_rdc = reduce_2_etope(E)
% E_rdc = reduce_2_etope(E,n_gen_to_remove,n_con_to_keep)
% E_rdc = reduce_2_etope(E,n_gen_to_remove,n_con_to_keep,flag_force_reduce)
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
% Updated: 16 Mar 2022

%% setup
    % sanity check
    if E.p_norm ~= 2
        error('This function only works for 2-ellipsotopes!')
    end
    
    % set default number of generators to reduce
    if nargin < 2
        n_gen_to_remove = 1 ;
    end
    
    if nargin < 3
        n_con_to_keep = 3 ; % threeee is a magic number
    end
    
    if nargin < 4
        flag_force_reduce = false ;
    end
    
    % desired number of gennies
    n_gen_orig = E.n_generators ;
    n_des = n_gen_orig - n_gen_to_remove ;
    
%% check if basic, then things are easy-peasy
    if E.is_basic()
        G_rdc = reduce_2_etope_generator_matrix(E.generators) ;
        E = ellipsotope(E.p_norm,E.center,G_rdc) ;
        return ;
    end
    
%% get into minimal form
    E = reduce_2_etope_to_minimal_exact_rep(E) ;
    n_gen = E.n_generators ;
    
% %% try using component ellipsotopes
%     if n_gen > (n_gen_orig - n_gen_to_remove)
%         E = reduce_component_2_etopes(E,n_gen_to_remove) ;
%         
%         n_gen = E.n_generators ;
%     end
%     
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
    else
        return ;
    end
    
    n_gen = E.n_generators ;
    
%% try component (constrained) zonotopes
    if n_gen > n_des
        E = reduce_component_zonotope(E,n_gen - n_des) ;
    else
        return ;
    end
    
    n_gen = E.n_generators ;
    
%% force reduction
    if flag_force_reduce
        while (n_con > 0) && (n_gen > n_des)
            % we'll try both methods, then see which one worked better
            
            % reduce a constraint
            E_con = reduce_etope_constraint_and_generator(E) ;
            E_con = reduce_2_etope_to_minimal_exact_rep(E_con) ;
            
            % pop and reduce
            E_pop = pop_generator(E) ;
            E_pop = reduce_2_etope_to_minimal_exact_rep(E_pop) ;
            E_pop = reduce_component_zonotope(E_pop) ;
            
            if E_con.n_generators > E_pop.n_generators
                E = E_pop ;
            else
                E = E_con ;
            end
            
            n_gen = E.n_generators ;
            n_con = E.n_constraints ;
        end
        
        % the last resort is just pop-n-drop
        while (n_gen > n_des)
            E = pop_generator(E) ;
            E = reduce_2_etope_to_minimal_exact_rep(E) ;
            E = reduce_component_zonotope(E) ;
            n_gen = E.n_generators ;
        end
    end

%% final check
    n_gen = E.n_generators ;
    if n_gen > n_des
        warning(['The ellipsotope was not reduced as much as you wanted! ',newline,...
            'The final ellipsotope has ',num2str(n_gen),' generators, ',...
            'whereas you wanted ',num2str(n_des),'.',newline,...
            'Consider trying again and forcing order reduction by calling: ',...
            newline,newline,...
            '    reduce_2_etope(E,',num2str(n_gen_to_remove),',',num2str(n_con_to_keep),',true)',...
            newline])
    end
end
function [E_rdc,idx_rdc] = reduce_etope_constraint_and_generator(E,idx_rdc)
% E_rdc = reduce_etope_constraint_and_generator(E)
% E_rdc = reduce_etope_constraint_and_generator(E,idx_rdc
%
% This script implements the method in [1] for eliminating one constraint
% and one generator of a constrained zonotope, which also works for
% 2-ellipsotopes.
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% See also: estimate_interval_ball_product_intersection.m,
%           estimate_best_idx_for_gen_and_con_removal.m
%
% Authors: Shreyas Kousik
% Created: 20 July 2021
% Updated: 16 March 2022 (added optional second input arg)

% get best index to reduce
if nargin < 2
    idx_rdc = estimate_best_idx_for_gen_and_con_removal(E) ;
end

% reduce
E_rdc = reduce_etope_constraint_and_generator_helper(E,idx_rdc) ;

% cleanup!
E_rdc.clean_properties() ;
end

%% helper function
function E_rdc = reduce_etope_constraint_and_generator_helper(E,idx_rdc)

    [p_norm,c,G,A,b,I,~,n_gen,n_con] = E.get_properties() ;

    % find a constraint where the j-th entry is nonzero (this must always
    % exist or else we can't reduce that constraint! ha ha ha ha)
    a_ij_all = A(:,idx_rdc) ;
    idx_con_rdc = find(a_ij_all ~= 0,1) ;
    
    % solve the first constraint where a_ij is nonzero for the coef value
    E_ij = zeros(n_gen,n_con) ;
    E_ij(idx_rdc,idx_con_rdc) = 1 ;
    a_ij_inv = 1./(A(idx_con_rdc,idx_rdc)) ;

    % create Gamma and Lambda as in Prop. 5 of [1]
    Gm = G*E_ij*a_ij_inv ;
    Lm = A*E_ij*a_ij_inv ;

    % create new etope parameters
    c_rdc = c + Gm*b ;
    G_rdc = G - Gm*A ;
    A_rdc = A - Lm*A ;
    b_rdc = b - Lm*b ;

    % delete i-th constraint and j-th generator column
    A_rdc(idx_con_rdc,:) = [] ;
    b_rdc(idx_con_rdc) = [] ;
    A_rdc(:,idx_rdc) = [] ;

    % delete j-th generator
    G_rdc(:,idx_rdc) = [] ;

    % reorganize the index set
    n_I = length(I) ;
    I_rdc = cell(1,n_I) ;
    for idx = 1:n_I
        J = I{idx} ;
        J_log = J == idx_rdc ;

        if any(J_log)
            J(J_log) = [] ;
        end

        J_log_minus = J > idx_rdc ;
        J(J_log_minus) = J(J_log_minus) - 1 ;

        I_rdc{idx} = J ;
    end
    
    % delete any empty index subsets
    I_empty_log = cellfun(@isempty,I_rdc) ;
    I_rdc = I_rdc(~I_empty_log) ;
    
    % clean up any empty constraints
    if isempty(A_rdc)
        A_rdc = [] ;
        b_rdc = [] ;
    end

    E_rdc = ellipsotope(p_norm,c_rdc,G_rdc,A_rdc,b_rdc,I_rdc) ;
end
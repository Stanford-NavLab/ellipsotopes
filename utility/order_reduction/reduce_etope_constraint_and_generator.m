function E_rdc = reduce_etope_constraint_and_generator(E)
% E_rdc = reduce_etope_constraint_and_generator(E)
%
% This script implements the method in [1] for eliminating one constraint
% and one generator of a constrained zonotope, which also works for
% 2-ellipsotopes.
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% See also: estimate_interval_ball_product_intersection.m
%
% Authors: Shreyas Kousik
% Created: 20 July 2021
% Updated: 14 March 2022 (functionized some stuff and fixed some bugz)

% get properties
%[p,~,G,A,b,~,~,n_gen,n_con,~] = get_properties(E) ;
[p,c,G,A,b,I,n_dim,n_gen,n_con,n_I] = get_properties(E) ;

if p ~= 2
    error('This function only works for 2-ellipsotopes for now!')
end

if isempty(A)
    error('The ellipsotope is unconstrained!')
end

%% compute bounds on ball product intersecting affine subspace
[~,~,int_R_rescaled] = estimate_interval_ball_product_intersection(E) ;

%% test computing approximate Hausdorff distance
% this applies the method in the appendix of [1]

% preallocate for solutions
H_hat = nan(1,n_gen) ;

% compute r_j values as in (A.6) of [1]
r_j_all = max([zeros(n_gen,1), max(abs(int_R_rescaled),[],2) - 1],[],2) ;

% check if any r_j = 0 (within some numerical tolerance)
r_j_test = abs(r_j_all) < 1e-8 ;
if any(r_j_test)
    idx_j_min = find(r_j_test,1) ;
else
    % compute H_hat approx for each j
    for idx_j = 1:n_gen
        % create j unit vector
        e_j = zeros(n_gen,1) ;
        e_j(idx_j) = 1 ;

        % create QP solution matrix as in (A.7)
        LHS = [G'*G + eye(n_gen), A', e_j ;
            A, zeros(n_con,n_con), zeros(n_con,1) ;
            e_j', zeros(1,n_con), 0] ;
        RHS = [zeros(n_gen+n_con,1) ; r_j_all(idx_j)] ;

        % solve for d_hat and lambda_hat (Shreyas was lazy and didn't
        % pre-factorize the Q = (G'G + I) matrix, but like, this could be done
        % faster in the future I guess...
        sol = pinv(LHS)*RHS ;

        % extract d_hat
        d_hat = sol(1:n_gen) ;

        % compute H_hat approximate cost
        H_hat(idx_j) = vecnorm(G*d_hat,2) + vecnorm(d_hat) ;
    end

    [~,idx_j_min] = min(H_hat) ;
end

%% reduce
E_rdc = reduce_one_con_and_gen(E,idx_j_min) ;

% cleanup!
E_rdc.clean_properties() ;

end

%% helper functions
function E_rdc = reduce_one_con_and_gen(E,j_rdc)
    [p_norm,c,G,A,b,I,~,n_gen,n_con] = E.get_properties() ;

    % find a constraint where the j-th entry is nonzero (this must always
    % exist or else we can't reduce that constraint! ha ha ha ha)
    a_ij_all = A(:,j_rdc) ;
    i_rdc = find(a_ij_all ~= 0,1) ;
    
    % solve the first constraint where a_ij is nonzero for the coef value
    E_ij = zeros(n_gen,n_con) ;
    E_ij(j_rdc,i_rdc) = 1 ;
    a_ij_inv = 1./(A(i_rdc,j_rdc)) ;

    % create Gamma and Lambda as in Prop. 5 of [1]
    Gm = G*E_ij*a_ij_inv ;
    Lm = A*E_ij*a_ij_inv ;

    % create new etope parameters
    c_rdc = c + Gm*b ;
    G_rdc = G - Gm*A ;
    A_rdc = A - Lm*A ;
    b_rdc = b - Lm*b ;

    % delete i-th constraint and j-th generator column
    A_rdc(i_rdc,:) = [] ;
    b_rdc(i_rdc) = [] ;
    A_rdc(:,j_rdc) = [] ;

    % delete j-th generator
    G_rdc(:,j_rdc) = [] ;

    % reorganize the index set
    n_I = length(I) ;
    I_rdc = cell(1,n_I) ;
    for idx = 1:n_I
        J = I{idx} ;
        J_log = J == j_rdc ;

        if any(J_log)
            J(J_log) = [] ;
        end

        J_log_minus = J > j_rdc ;
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
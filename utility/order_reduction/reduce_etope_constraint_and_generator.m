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
% Authors: Shreyas Kousik
% Created: 20 July 2021
% Updated: --

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
% initialize
int_E = [-ones(n_gen,1), ones(n_gen,1)] ;
int_R = [-inf(n_gen,1), inf(n_gen,1)] ;

% iterate Algorithm 1 of [1]
for idx_i = 1:n_con
    for idx_j = 1:n_gen
        a_ij = A(idx_i,idx_j) ;
        
        if a_ij ~= 0
            R_j = int_R(idx_j,:) ;
            E_j = int_E(idx_j,:) ;
            
            a_ij_inv = 1./a_ij ;
            
            idxs_k = 1:n_gen ;
            idxs_k(idx_j) = [] ;
                        
            R_j_RHS = a_ij_inv*b(idx_i).*[1 1] ;
            E_k_sum = [0 0] ;
            for idx_k = idxs_k
                a_ik = A(idx_i,idx_k) ;
                E_k = int_E(idx_k,:) ;
                E_k_prod = interval_scalar_mult(a_ij_inv*a_ik,E_k) ;
                E_k_sum = interval_add(E_k_sum,E_k_prod) ;
            end
            
            R_j_RHS = interval_subtract(R_j_RHS, E_k_sum) ;
            
            R_j = interval_intersect(R_j, R_j_RHS) ;
            
            E_j = interval_intersect(E_j,R_j) ;
            
            int_R(idx_j,:) = R_j ;
            int_E(idx_j,:) = E_j ;
        end
    end
end

%% test computing approximate Hausdorff distance
% this applies the method in the appendix of [1]

% preallocate for solutions
H_hat = nan(1,n_gen) ;

% compute r_j values
int_r = max(abs(int_R),[],2) ;

% compute H_hat approx for each j
for idx_j = 1:n_gen
    % create j unit vector
    e_j = zeros(n_gen,1) ;
    e_j(idx_j) = 1 ;
       
    % create matrix
    LHS = [G'*G + eye(n_gen), A', e_j ;
        A, zeros(n_con,n_con), zeros(n_con,1) ;
        e_j', zeros(1,n_con), 0] ;
    RHS = [zeros(n_gen+n_con,1) ; int_r(idx_j)] ;
    
    % solve for d_hat and lambda_hat
    sol = inv(LHS)*RHS ;
    
    % extract d_hat
    d_hat = sol(1:n_gen) ;
    
    % compute H_hat approximate cost
    H_hat(idx_j) = vecnorm(G*d_hat,2) + vecnorm(d_hat) ;
end

[~,idx_j_min] = min(H_hat) ;

%% reduce
E_rdc = reduce_one_con_and_gen(E,idx_j_min) ;

% cleanup!
E_rdc.clean_properties() ;

end

%% helper functions
function I = interval_intersect(I_1,I_2)
     I_lo = [I_1(:,1), I_2(:,1)] ;
     I_hi = [I_1(:,2), I_2(:,2)] ;
     
     I = [max(I_lo,[],2), min(I_hi,[],2)] ;
end

function I = interval_add(I_1,I_2)
    I = [I_1(:,1) + I_2(:,1), I_1(:,2) + I_2(:,2)] ;
end

function I = interval_subtract(I_1,I_2)
    I = [I_1(:,1) - I_2(:,2), I_1(:,2) - I_2(:,1)] ;
end

function I = interval_scalar_mult(s,I)
    I = s.*I ;
    if s < 0
        I = [I(:,2), I(:,1)] ;
    end
end

function E_rdc = reduce_one_con_and_gen(E,j_rdc)
    [p_norm,c,G,A,b,I,~,n_gen,n_con] = E.get_properties() ;

    % solve the first constraint for the coefficient value
    E_j1 = zeros(n_gen,n_con) ;
    E_j1(j_rdc,1) = 1 ;
    a_1j_inv = 1./A(1,j_rdc) ;

    % create Gamma and Lambda
    Gm = G*E_j1*a_1j_inv ;
    Lm = A*E_j1*a_1j_inv ;

    % create new etope
    c_rdc = c + Gm*b ;
    G_rdc = G - Gm*A ;
    A_rdc = A - Lm*A ;
    b_rdc = b - Lm*b ;

    % delete first constraint and j-th generator column
    A_rdc = A_rdc(2:end,:) ;
    b_rdc = b_rdc(2:end) ;
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
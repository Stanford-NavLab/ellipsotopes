function [int_E,int_R,int_R_rescaled] = estimate_interval_ball_product_intersection(E)
% [int_E,int_R] = estimate_interval_ball_product_intersection(E)
%
% This script implements the method in Alg. 1 of [1] for computing an
% interval containing the intersection of the the B_infty norm ball and the
% hyperplane defined by A*xi = b.
%
% The output int_E is an interval subset of the B_infty ball
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% Authors: Shreyas Kousik
% Created: 14 Mar 2022
% Updated: --

    [~,~,~,A,b,~,~,n_gen,n_con] = get_properties(E) ;
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
    
    % get interval scale
    xi_m = 0.5*(int_E(:,1) + int_E(:,2)) ;
    xi_r = 0.5*(int_E(:,2) - int_E(:,1)) ;
    
    % rescale R
    int_R_rescaled = (int_R - repmat(xi_m,1,2))./repmat(xi_r,1,2) ;
end

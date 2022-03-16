function idx_rdc = estimate_best_idx_for_gen_and_con_removal(E)
% idx_rdc = estimate_best_idx_for_gen_and_con_removal(E)
%
% Use the method in [1] for eliminating one constraint and one generator of
% an ellipsotope; this function identifies the index idx_rdc for use in
% order reduction, where idx_rdc is in 1,...,n_gen.
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% See also: reduce_etope_constraint_and_generator.m,
%           reduce_etope_constraint.m,
%           estimate_interval_ball_product_intersection.m
%
% Authors: Shreyas Kousik
% Created: 16 Mar 2022
% Updated: --
    [~,~,G,A,~,~,~,n_gen,n_con,~] = get_properties(E) ;

    if isempty(A)
        error('The ellipsotope is unconstrained!')
    end

    %% compute bounds on ball product intersecting affine subspace
    [~,~,int_R_rescaled] = estimate_interval_ball_product_intersection(E) ;

    %% test computing approximate Hausdorff distance
    if nargin < 2
        % this applies the method in the appendix of [1]

        % preallocate for solutions
        H_hat = nan(1,n_gen) ;

        % compute r_j values as in (A.6) of [1]
        r_j_all = max([zeros(n_gen,1), max(abs(int_R_rescaled),[],2) - 1],[],2) ;

        % check if any r_j = 0 (within some numerical tolerance)
        r_j_test = abs(r_j_all) < 1e-8 ;
        if any(r_j_test)
            idx_rdc = find(r_j_test,1) ;
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

            [~,idx_rdc] = min(H_hat) ;
        end
    end
end
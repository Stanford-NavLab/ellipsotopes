function [G_MVOE,Q_MVOE] = make_MVOE_generator_matrix(G_1,G_2)
% [G_MVOE,Q_MVOE] = make_MVOE_generator_matrix(G_1,G_2)
%
% Compute the shape matrix for the minimum volume outer ellipsoid of the
% Minkowski sum of the pair of ellipsoids defined by Q_1 and Q_2 using the
% technique in [1].
%
% [1] Halder, A., 2018, December. On the parameterized computation of
%     minimum volume outer ellipsoid of Minkowski sum of ellipsoids. In
%     2018 IEEE Conference on Decision and Control (CDC) (pp. 4040-4045).
%     IEEE.
%
% Authors: Shreyas Kousik
% Created: 13 June 2021
% Updated: --
    
    % tolerate for iteration
    tol = 1e-10 ;

    % get shape matrices
    Q_1 = inv(pinv(G_1)'*pinv(G_1)) ;
    Q_2 = inv(pinv(G_2)'*pinv(G_2)) ;

    % get R matrix
    R = inv(Q_1)*Q_2 ;

    % get lambdas (eigenvalues)
    lm = eig(R) ;

    % perform iteration until convergence
    bt = 0 ;
    tol_violation = varphi(bt,lm) ;

    while tol_violation > tol
        % iterate eq. (21) of [1]
        n = sum(1./(1 + bt.*lm)) ;
        d = sum(lm./(1 + bt.*lm)) ;
        bt = sqrt(n/d) ;
        tol_violation = varphi(bt,lm) ;
    end

    % construct new ellipsotope shape matrix
    Q_MVOE = inv((1 + 1/bt).*Q_1 + (1 + bt).*Q_2) ;
    G_MVOE = inv(Q_MVOE^(1/2)) ;
end

%% helper functions
function val = varphi(b,lm)
   val = sum((1 - (b.^2).*lm)./(1 + b.*lm)) ;
end
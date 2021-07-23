function out = test_point_in_ball_product(p_norm,P,I)
% out = test_point_in_ball_product(p_norm,P,I)
%
% Test is the point/s P live in the ball product defined by p_norm and I.
%
% Authors: Shreyas Kousik
% Created: 23 July 2021
    n_I = length(I) ;
    n_P = size(P,2) ;

    out = false(n_I,n_P) ;

    for idx = 1:n_I
        % get index subset
        J = I{idx} ;

        % test if points obey p-norm
        V = vecnorm(P(J,:),p_norm,1) ;
        out(idx,:) =  all(V <= 1,1) ;
    end

    out = all(out,1) ;
end
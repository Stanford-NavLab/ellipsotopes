function P = sample_from_ellipsotope(E,n_P)
% P = sample_from_ellipsotope(E)
% P = sample_from_ellipsotope(E,n_points)
%
% Sample random point/s from inside the ellipsotope E. By default this
% outputs a single point.
%
% Authors: Shreyas Kousik
% Created: 23 July 2021

    if nargin < 2
        n_P = 1 ;
    end

    % get properties
    [p_norm,c,G,A,b,I,~,n_gen] = get_properties(E) ;

    % set up how many points to sample initially
    n_B = max(10^5,n_P) ;

    % create a bunch of points that are feasible to the linear constraints if
    % they exist
    if ~isempty(A)
        % get a feasible point
        t = pinv(A)*b ;

        % get the nullspace of the constraint
        K = null(A) ;
        n_K = size(K,2) ;

        % create a bunch of points that satisfy the constraint
        d_K = 4.*rand(n_K,n_B) - 2 ;
        B = K*d_K + repmat(t,1,n_B) ;
    else
        % sample from infinity norm ball
        B = 2.*rand(n_gen,n_B) - 1 ;
    end

    % test which of the points satisfy the p-norm ball constraint
    B_log = test_point_in_ball_product(p_norm,B,I) ;

    % keep the points that were ok
    B = B(:,B_log) ;
    n_B = size(B,2) ;

    if n_B > n_P
        B = B(:,1:n_P) ;
    else
        n_P = n_B ;
    end

    % map points to workspace
    P = repmat(c,1,n_P) + G*B ;

end
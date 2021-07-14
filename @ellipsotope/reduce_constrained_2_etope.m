function E_rdc = reduce_constrained_2_etope(E)

% get propz
[p,c,G,A,b,~,~,~,~,n_I] = get_properties(E) ;

if n_I > 1
    error('Cannot use this method for an indexed 2-etope!')
else
    % get affine map to hyperplane
    t = pinv(A)*b ; % center of intersection area
    rd = sqrt(1 - vecnorm(t)^2) ; % radius of intersected ball
    K = null(A) ;
    T = rd*K ;

    % compute new etope matrix
    c_rdc = c + G*t ;
    G_rdc = G*T ;

    % reduce G_rdc to 2-D
    G_rdc = reduce_2_etope_generator_matrix(G_rdc) ;

    % create new tope
    E_rdc = ellipsotope(2,c_rdc,G_rdc) ;
end
end
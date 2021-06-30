function E_rot = rotated_body_ellip(robot)
    
    G = robot.rot_body_zono.generators; 
    n_gen = size(G,2); n_dim = size(G,1);
    % compute E_0^(-1/2)
    E_0 = n_gen*(G*G') ;
    E_0_inv_sqrt = E_0^(-1/2) ;
    G_t = E_0_inv_sqrt*G ;
    % solve [1,(19)] with YALMIP
    lm = sdpvar(n_gen,1) ; % decision variable
    options = sdpsettings('verbose',0) ; 
    obj = sum(lm) ;
    cons = [diag(lm) - G_t'*G_t >= 0, lm >= 0] ;
    optimize(cons,obj,options) ;
    r_hat = sum(value(lm)) ;
    % compute enclosing ellipsoid
    E_enc = r_hat * E_0 ;
    E_rot = ellipsoid(E_enc,zeros(n_dim,1)) ;

end